# ---------------------------- DESCRIPTION --------------------------------------- #

"""
Function to compute the master flat field from a set of flat-field observations of 
a single observation mode. 

Modified by FJBM on 2025-10-28.
2nd version (2026-02-24): changed zero blueshift to coincide with minimum blueshift across FOV
Patched by JCT on 2026-04-15:
- support N_mods != 4
- map 0s -> 1 and 3.02 -> 2.02 for prefilter and blueshift guess

Modified by FJBM on 2026-06-19.
- Support all observation modes for blueshift and prefilter removal.
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Standard library
import time
from pathlib import Path
import logging
import os
from matplotlib import pyplot as plt

# Third-party libraries
import numpy as np
from tqdm import tqdm
from scipy.interpolate import interp1d,CubicSpline, PchipInterpolator
from scipy.optimize import minimize

# Local application modules
import config as cf
import image_handler as ih
import prefilter_removal_fjbm as pr
from joblib import Parallel, delayed

# ------------------------------  CODE  ------------------------------------------ # 

def compute_master_flat_field(flat_fields_paths, dc, lambda_repeat=4, verbose=False,
                              norm_method="blueshift", remove_prefilter=True,
                              pref_model=True, import_blueshift_guess=True,
                              modify_linearity=([1539, 1540], [1.0, 1.0]),
                              norm_roi=[300, -300, 300, -300], corte = 650, centro = [800,800],
                              discard_repetitions = False, flat_optimization = False):

    """
    Function to compute the flat-field observation from the images paths. 

    inputs:
        - flat_field_paths (list) : List containing all the paths to the images composing a single 
          flat-field observation. 
        - dc (np.array) : Dark current with shape (2, nx, ny)
        - lambda_repeat (int, default : 4) : Lambda repeat parameter of the observation.
        - norm_method (str, default : "blueshift") : Normalization method. avg, mod or blueshift.
        - remove_prefilter (bool, default : True) : Option to remove prefilter from the flats profiles. 
        - pref_model : If True uses analytic model in prefilter_removal_fjbm.py;
                       otherwise fits prefilter from flat profile.
        - volts (unused here / kept for compatibility)
        - import_blueshift_guess (bool, default : True) : Option to import previous
          blueshift guess to speed up the fitting process.
    returns:
        - ff_norm (np.array) : Array containing the flat field (cams, Nlambda, Nmods, Nx, Ny)
        - ff_info (dictionary) : Dictionary containing all info of the flat-field.
    """
    
    tic = time.time()
    block_size = 50
    print("CPUs:", os.cpu_count())

    if verbose:
        print(f"\nComputing flats.")
        print(f"------------------")
        print(f"N flats: {len(flat_fields_paths)}")

    # Read first image to get info from header
    _, h = ih.read(flat_fields_paths[0])
    om = h["ObservationMode"]

    # Get observation mode parameters
    N_wls = cf.om_config[om]["Nlambda"]
    N_mods = cf.om_config[om]["Nmods"]
#    wvlv = np.array(cf.om_config[om]["lambda_array"]) * 1e-13  # Wavelengths in meters

    # Calculate number of repetitions
    if len(flat_fields_paths) % (2 * N_wls * N_mods * lambda_repeat) == 0:
        nreps = int(len(flat_fields_paths) / (2 * N_wls * N_mods * lambda_repeat))
    else:
        raise Exception("Observations are incomplete, please remove images from incomplete OC. This will be upgraded...")
    print(f"Available repetitions: {nreps}")

    if discard_repetitions is None or discard_repetitions is False:
        discard_repetitions = []
    elif np.isscalar(discard_repetitions):
        discard_repetitions = [int(discard_repetitions)]
    else:
        discard_repetitions = list(discard_repetitions)

    for rep in discard_repetitions:
        if rep < 0 or rep >= nreps:
            raise ValueError(
                f"discard_repetitions contains invalid repetition {rep}. "
                f"Valid range is [0, {nreps-1}]"
            )

    valid_reps = [r for r in range(nreps) if r not in discard_repetitions]

    if verbose:
        print(f"Discarding repetitions: {discard_repetitions}")
        print(f"Using repetitions: {valid_reps}")

    n_files_rep = 2 * N_wls * N_mods * lambda_repeat

    flat_fields_paths = [
        path
        for rep in range(nreps)
        if rep not in discard_repetitions
        for path in flat_fields_paths[
            rep*n_files_rep:(rep+1)*n_files_rep
        ]
    ]

    nreps = len(valid_reps)

    if verbose:
        print(f"Observation Mode: {om}")
        print(f"Nº of repetitions: {nreps}")
        print(f"Nº of wavelengths: {N_wls}")
        print(f"Nº of Modulations: {N_mods}")
#        print(f"Wavelength sampling (pm): {wvlv * 1e12}")

    # -------------------------------------------------------------------------
    # Read images and dark-correct
    # -------------------------------------------------------------------------
    flat_obs = ih.nominal_flat(om, flat_fields_paths, nreps, dc, modify_linearity=modify_linearity)
    data = flat_obs.get_data()

    # -------------------------------------------------------------------------
    # Read commanded voltages and reorder if decreasing
    # -------------------------------------------------------------------------
    volts_list = [] 
    for lambd in range(N_wls):
            volts_list.append(flat_obs.info["Images_headers"][f"wv_{lambd}"][f"Mod_{0}"]['hvps_comm_volts'][0])

    volts_list = np.asarray(volts_list)
    dv = np.diff(volts_list)
    reordering = False
    if np.all(dv > 0):
        logging.info("Voltages increasing")
    elif np.all(dv < 0):
        logging.info("Voltages decreasing -> reordering")
        idx = np.argsort(volts_list)
        volts_list = volts_list[idx]
        data = data[:, idx, ...]
        reordering = True
    else:
        logging.warning("Voltages are not monotonic")

    # -------------------------------------------------------------------------
    # Convert volts to wavelength
    # -------------------------------------------------------------------------
    wvl=pr.volts_2_lambda(volts_list, pr.Config[cf.om_config[om]["line"]])*1e-10 #Wavelength in meters
    wvlv = wvl - pr.wvl0(om)  # Wavelengths in meters relative to central wavelength

    # -------------------------------------------------------------------------
    # Define cropped area for the blueshift fit
    # -------------------------------------------------------------------------
    x0 = centro[0] - corte
    y0 = centro[1] - corte
    xf = centro[0] + corte
    yf = centro[1] + corte
    Nx = xf - x0
    Ny = yf - y0
    logging.info(f'  >> flat (processed area) cropping..........{x0,y0,xf,yf,Nx,Ny}')

    # Region where blueshift is smallest (ad-hoc)
    x0_mean = 1200
    y0_mean = 1200
    xf_mean = x0_mean + 100
    yf_mean = y0_mean + 100

    # -------------------------------------------------------------------------
    # Average profiles in the reference region
    # -------------------------------------------------------------------------
    cam1_ave = np.mean(data[0, :, :, x0_mean:xf_mean, y0_mean:yf_mean], axis=(1, 2, 3))
    cam2_ave = np.mean(data[1, :, :, x0_mean:xf_mean, y0_mean:yf_mean], axis=(1, 2, 3))

    norm_ff = np.ones_like(data)

    # -------------------------------------------------------------------------
    # Normalize flat field according to selected method
    # -------------------------------------------------------------------------
    if norm_method == "blueshift":
        # Interpolation of the averaged profiles over wavelength TODO: se puede usar otro no?
        cam1_interp_ref = interp1d(wvlv, cam1_ave, kind='quadratic', bounds_error=False,
                               fill_value='extrapolate')
        cam2_interp_ref = interp1d(wvlv, cam2_ave, kind='quadratic', bounds_error=False,
                               fill_value='extrapolate')

        meth = 'Nelder-Mead'

        # ---------------------------------------------------------------------
        # Blueshift initial guess
        # ---------------------------------------------------------------------
        if import_blueshift_guess:
            if om == "0s" or om == "0p" or om=="1" or om=="4":
                om_guess = "1"
            elif om== "2.02" or om == "3.02" or om=="5.02":
                om_guess = "2.02"
            elif om == "2.06" or om == "3.06" or om=="5.06":
                om_guess = "2.02" #For 525.06 we use the same guess as for 525.02

            if om_guess == '2.02':
                fname_guess = 'ff_om2.02_D10-4340-5619_nonorm'
            elif om_guess == '1':
                fname_guess = 'ff_om1_D10-2740-4339_nonorm'
            else:
                raise ValueError(f"No blueshift guess configured for ObservationMode={om} (mapped to {om_guess})")

            if verbose:
                logging.info(f"Using blueshift guess mode: {om_guess} for ObservationMode {om}")

            BASE_DIR = Path(__file__).resolve().parent
            blueshift_guess = np.load(BASE_DIR / f"{fname_guess}_wvl_shifts_cam1.npy")
            scale_guess = np.load(BASE_DIR / f"{fname_guess}_scales_cam1.npy")

        # Output arrays

        wvl_shifts = np.zeros((2, Nx, Ny))
        scales = np.zeros((2, Nx, Ny))

        logging.info('Fitting blueshift and scale for each camera...')

        for cam in range(2):
            logging.info(f'Camera {cam + 1}')
            cam_ave = cam1_ave if cam == 0 else cam2_ave

            interp_method = 'spline'   # 'quadratic', 'spline', 'pchip'
            plot_interpolators = True

            # ---------------------------------------------------
            # Preparar datos
            # ---------------------------------------------------

            if om == '2.02':
                wvlv_interp = np.append(wvlv, 30e-12)
                cam_ave_interp = np.append(cam_ave, cam_ave[-1])
            else:
                wvlv_interp = wvlv.copy()
                cam_ave_interp = cam_ave.copy()

            interp_sort_idx = np.argsort(wvlv_interp)
            x_interp = wvlv_interp[interp_sort_idx]
            y_interp = cam_ave_interp[interp_sort_idx]

            # ---------------------------------------------------
            # Crear interpolador seleccionado
            # ---------------------------------------------------
            if interp_method.lower() == 'quadratic':
                line_interp = interp1d(
                    x_interp, y_interp,
                    kind='quadratic',
                    bounds_error=False,
                    fill_value='extrapolate'
                )
            elif interp_method.lower() == 'pchip':
                line_interp = PchipInterpolator(
                    x_interp, y_interp, extrapolate=True
                )
            else:
                line_interp = CubicSpline(
                    x_interp, y_interp,
                    bc_type='natural',
                    extrapolate=True
                )

            # ---------------------------------------------------
            # Diagnóstico opcional
            # ---------------------------------------------------
            if plot_interpolators:
                wvlv2 = np.linspace(np.min(x_interp), np.max(x_interp), 5000)

                interp_quad = interp1d(
                    x_interp, y_interp, kind='quadratic',
                    bounds_error=False, fill_value='extrapolate'
                )
                interp_cubic = CubicSpline(
                    x_interp, y_interp,
                    bc_type='natural',
                    extrapolate=True
                )
                interp_pchip = PchipInterpolator(
                    x_interp, y_interp, extrapolate=True
                )

                prof_quad = interp_quad(wvlv2)
                prof_cubic = interp_cubic(wvlv2)
                prof_pchip = interp_pchip(wvlv2)

                cen_quad = wvlv2[np.nanargmin(prof_quad)]
                cen_cubic = wvlv2[np.nanargmin(prof_cubic)]
                cen_pchip = wvlv2[np.nanargmin(prof_pchip)]

                plt.figure(figsize=(9, 5))
                plt.plot(x_interp, y_interp, 'ko', label='Data')
                plt.plot(wvlv2, prof_quad, label=f'Quadratic ({cen_quad*1e12:.3f} pm)')
                plt.plot(wvlv2, prof_cubic, label=f'Spline ({cen_cubic*1e12:.3f} pm)')
                plt.plot(wvlv2, prof_pchip, label=f'PCHIP ({cen_pchip*1e12:.3f} pm)')
                plt.legend()
                plt.grid(True)
                plt.show()

                print(f"Quadratic : {cen_quad*1e12:.4f} pm")
                print(f"Spline    : {cen_cubic*1e12:.4f} pm")
                print(f"PCHIP     : {cen_pchip*1e12:.4f} pm")

            # ---------------------------------------------------
            # Centro de línea con el interpolador elegido
            # ---------------------------------------------------

            wvlv2 = np.linspace(np.min(x_interp), np.max(x_interp), 5000)

            line_center = wvlv2[np.nanargmin(line_interp(wvlv2))]

            ind_fit0 = 0
            ind_fitf = -1

            wvlv_fit = wvlv[ind_fit0:ind_fitf]
            reference_interp = cam1_interp_ref if cam == 0 else cam2_interp_ref

            # ---------------------------------------------
            # Fit blueshift and scale
            # ---------------------------------------------
            if flat_optimization:
                logging.info("using opmtimized version")
                
                # Rejilla de shifts
                delta_grid = np.linspace(-10e-12, 10e-12, 510) #increase

                # Referencias interpoladas precalculadas
                refs = np.array([
                    reference_interp(wvlv_fit + delta)
                    for delta in delta_grid
                ])
                refs_norm = np.sum(refs**2, axis=1)

                def fit_block(i0, i1):
                    shifts_block = np.empty((i1 - i0, Ny))
                    scales_block = np.empty((i1 - i0, Ny))
                    for ii, i in enumerate(range(i0, i1)):
                        for j in range(Ny):
                            cam_indiv = np.mean(
                                data[cam, ind_fit0:ind_fitf, :, x0+i, y0+j],
                                axis=1
                            )
                            dots = refs @ cam_indiv
                            scales_test = dots / refs_norm
                            scales_test = np.clip(scales_test, 0.1, 2.0)
                            residuals = refs * scales_test[:, None] - cam_indiv[None, :]
                            merits = np.sum(residuals**2, axis=1)
                            best = np.argmin(merits)
                            shifts_block[ii, j] = delta_grid[best]
                            scales_block[ii, j] = scales_test[best]
                    return i0, i1, shifts_block, scales_block

                # --------------------------------------------------------
                # Paralelización
                # --------------------------------------------------------
                blocks = [
                    (i0, min(i0 + block_size, Nx))
                    for i0 in range(0, Nx, block_size)
                ]
                results = Parallel(
                    n_jobs=-1,
                    backend="loky",
                    verbose=10
                )(
                    delayed(fit_block)(i0, i1)
                    for i0, i1 in blocks
                )
                # --------------------------------------------------------
                # Guardar resultados
                # --------------------------------------------------------
                for i0, i1, shifts_block, scales_block in results:

                    wvl_shifts[cam, i0:i1, :] = shifts_block
                    scales[cam, i0:i1, :] = scales_block

            else:
                # Fit each pixel
                for i in tqdm(range(Nx)):
                    for j in range(Ny):
                        cam_indiv = np.mean(data[cam, ind_fit0:ind_fitf, :, x0+i, y0+j], axis=1)

                        def merit_indiv1(params):
                            delta_wvl = params[0]
                            scale = params[1]

                            wvlv_shifted = wvlv_fit + delta_wvl
                            cam_mean_shifted = reference_interp(wvlv_shifted)

                            diff = scale * cam_mean_shifted - cam_indiv
                            return np.sum(diff ** 2)

                        if import_blueshift_guess:
                            guess_ij = [blueshift_guess[i, j], scale_guess[i, j]]
                        else:
                            guess_ij = [0, 1]

                        minim = minimize(
                            merit_indiv1,
                            x0=guess_ij,
                            method=meth,
                            bounds=[(-10e-12, 10e-12), (0.1, 2)]
                        )

                        wvl_shifts[cam, i, j] = minim.x[0]
                        scales[cam, i, j] = minim.x[1]

            # ---------------------------------------------
            # Re-reference blueshift
            # Blueshift in the Sun is negative, so using max_shift as reference
            # makes the least-shifted pixel equal to 0.
            # Example: [-5,-4,-3,-2] -> subtract(-2) -> [-3,-2,-1,0]
            # ---------------------------------------------
            raw_min_shift = np.min(wvl_shifts[cam])
            raw_max_shift = np.max(wvl_shifts[cam])

            logging.info(
                f"Camera {cam+1}: raw blueshift range = "
                f"[{raw_min_shift*1e12:.3f}, {raw_max_shift*1e12:.3f}] pm"
            )

            max_shift = raw_max_shift
            wvl_shifts[cam] -= max_shift

            logging.info(
                f"Camera {cam+1}: referenced blueshift range = "
                f"[{np.min(wvl_shifts[cam])*1e12:.3f}, {np.max(wvl_shifts[cam])*1e12:.3f}] pm"
            )

            wvl_offset = line_center + max_shift

            if verbose:
                logging.info(
                    f'Sampling offset with respect to line center for cam {cam + 1} (pm): {wvl_offset * 1e12:g}'
                )
            # ---------------------------------------------
            # Apply Eq. 3 / correct profiles
            # ---------------------------------------------
            def correct_block(i0, i1):
                logging.info(f"Starting block {i0}:{i1}")

                norm_block = np.empty((i1 - i0, N_wls, N_mods, Ny))

                for ii, i in enumerate(range(i0, i1)):
                    for j in range(Ny):
                        cam_shifted_average = np.zeros(N_wls)

                        for k in range(N_mods):
                            if om == '2.02':
                                cam_data_interp = np.append(
                                    data[cam, :, k, x0+i, y0+j],
                                    data[cam, -1, k, x0+i, y0+j]
                                )
                            else:
                                cam_data_interp = data[cam, :, k, x0+i, y0+j]

                            # Always sort interpolation axis
                            cam_data_interp_sorted = cam_data_interp[interp_sort_idx]

                            corr_interp = interp1d(
                                x_interp,
                                cam_data_interp_sorted,
                                kind='quadratic',
                                bounds_error=False,
                                fill_value='extrapolate'
                            )

                            cam_shifted = corr_interp(
                                wvlv - wvl_shifts[cam, i, j]
                            )
                            cam_shifted_average += cam_shifted

                        cam_shifted_average /= N_mods

                        for k in range(N_mods):
                            norm_block[ii, :, k, j] = (
                                scales[cam, i, j]
                                * data[cam, :, k, x0+i, y0+j]
                                / cam_shifted_average
                            )

                logging.info(f"Block {i0}:{i1} finished")
                return i0, i1, norm_block

            blocks = [
                (i0, min(i0 + block_size, Nx))
                for i0 in range(0, Nx, block_size)
            ]

            results = Parallel(
                n_jobs=-1,
                backend="loky",
                verbose=10
            )(
                delayed(correct_block)(i0, i1)
                for i0, i1 in blocks
            )

            for i0, i1, norm_block in results:
                norm_ff[cam, :, :, x0+i0:x0+i1, y0:y0+Ny] = np.moveaxis(norm_block, 0, 2)

        np.savez('blueshift.npz', wvl_shifts=wvl_shifts, scales=scales)

    elif norm_method == "avg":
        norma = np.mean(
            data[:, :, :, norm_roi[0]:norm_roi[1], norm_roi[2]:norm_roi[3]],
            axis=(2, 3, 4)
        )
        for lambd in range(N_wls):
            for mod in range(N_mods):
                norm_ff[0, lambd, mod] = data[0, lambd, mod] / norma[0, lambd]
                norm_ff[1, lambd, mod] = data[1, lambd, mod] / norma[1, lambd]
    
    elif norm_method == "mod":
        norma = np.zeros(np.shape(data[:, :, :]))
        for lambd in range(N_wls):
            for mod in range(N_mods):
                norma[:, lambd, mod] = np.mean(
                    data[:, lambd, mod, norm_roi[0]:norm_roi[1], norm_roi[2]:norm_roi[3]],
                    axis=(1, 2)
                )
                norm_ff[0, lambd, mod] = data[0, lambd, mod] / norma[0, lambd, mod]
                norm_ff[1, lambd, mod] = data[1, lambd, mod] / norma[1, lambd, mod]

    elif norm_method == "none" or norm_method is None:
        norm_ff = data

    else:
        raise Exception("Invalid normalization method. Please select 'blueshift','avg', 'mod' or 'none'")

    debug = False
    if debug:
        from joblib import dump
        dump(locals(), "tmp/mydebug.joblib")
        raise SystemExit

    # -------------------------------------------------------------------------
    # Remove prefilter if required
    # -------------------------------------------------------------------------
    if remove_prefilter:
        om_pref = om

        if verbose:
            logging.info(f"Using prefilter mode: {om_pref} for ObservationMode {om}")

        if pref_model is True:
            prefilter = pr.prefilter_model(om_pref, wvlv)
        else:
            prefilter = pr.prefilter_fitting(cam1_ave, om_pref, wvlv)

        norm_ff *= prefilter[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis]

    if reordering:
        norm_ff = norm_ff[:, np.argsort(idx), ...]

    logging.info(f"Flat-fields computed in {round(time.time() - tic, 3)} s.")
    return norm_ff, flat_obs.get_info()


def correct_observation(data, ff, onelambda=False):
    """
    Function to apply the flat_field correction. 

    Inputs:
        - data (np.array) : observing mode data
        - ff (np.array) : flat_field data
        - onelambda (Boolean, default : False) : Select if only one lambda is passed.
    returns:
        - corrected : Corrected data.  
    """

    if onelambda:
        data = data[:, np.newaxis]

    shape = np.shape(data)
    nlambda = shape[1]
    nmods = shape[2]

    om_corr = np.zeros(np.shape(data))
    for lambd in range(nlambda):
        for mod in range(nmods):
            for cam in range(2):
                om_corr[cam, lambd, mod] = data[cam, lambd, mod] / ff[cam, lambd, mod]

    if onelambda:
        return om_corr[:, 0]
    else:
        return om_corr

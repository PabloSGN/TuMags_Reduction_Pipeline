# ---------------------------- DESCRIPTION --------------------------------------- #

"""
Function to compute the master flat field from a set of flat-field observations of 
a single observation mode. 
"""

# ------------------------------ IMPORTS ----------------------------------------- #
import numpy as np
import time

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

import scipy
from packaging.version import parse
# Import sims or simpson depending on scipy's version
if parse(scipy.__version__) >= parse("1.7"):  # Adjust the version as needed
    from scipy.integrate import simpson as simps
else:
    from scipy.integrate import simps

from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
from scipy.interpolate import RectBivariateSpline




# Own modules
import config as cf
import image_handler as ih
import prefilter_removal as pr

# ------------------------------  CODE  ------------------------------------------ # 

def compute_master_flat_field(flat_fields_paths, dc, lambda_repeat = 4, verbose = False, 
                              norm_method = "avg", remove_prefilter = False, pref_model = None, 
                              volts = None, compute_blueshift_flag = False, cpos = -1, plot_flag = False,
                              plotname = "blueshift.png"):
    """
    Function to compute the flat-field observation from the images paths. 

    inputs:
        - flat_field_paths (list) : List containing all the paths to the images composing a single 
        flat-field observation. 
        - dc (np.array) : Dark current. 
        - lambda_repeat (int, default : 4) : Lambda repeat parameter of the observation.
        - norm_method (str, default : "avg) : Normalization method. avg or mod.
        - remove_prefilter (Boolean, default : False) : Option to remove prefilter from the flats profiles. 
        - pref_model : Prefilter model file rerquired if remove_prefilter = True. 
        - volts (str / None, default = None) : Set to "read" if read voltages are to be used for the pref_removal.
        If None, fixed voltages are used. 
    returns:
        - ff_data (np.array) : Array containing the flat field (cams, Nlambda, Nmods, Nx, Ny)
        - ff_info (dictionary) : Dictionary containing all info of the flat-field.    
    """
    
    tic = time.time()

    if verbose:
        print(f"\nComputing flats.")
        print(f"------------------")
        print(f"N flats: {len(flat_fields_paths)}")

    # Read first image to get info from header
    _, h = ih.read(flat_fields_paths[0])

    om = h["ObservationMode"]
    line = cf.om_config[om]["line"]
    N_wls = cf.om_config[om]["Nlambda"]
    N_mods = cf.om_config[om]["Nmods"]

    if len(flat_fields_paths) % (2 * N_wls * N_mods * lambda_repeat) == 0:
        nreps = int(len(flat_fields_paths) / (2 * N_wls * N_mods * lambda_repeat))
    else:
        raise Exception("Observations are incomplete, please remove images from incomplete OC. This will be upgraded...")
    
    if verbose:
        print(f"Observation Mode: {om}")
        print(f"Nº of repetitions: {nreps}")
        print(f"Nº of wavelengths: {N_wls}")
        print(f"Nº of Modulations: {N_mods}")

    # Read images and correct them from dark current.
    flat_obs = ih.nominal_flat(om, flat_fields_paths, nreps, dc)

    data = flat_obs.get_data()

    # Normalize flat-fields
    norm_ff = np.zeros(np.shape(data))
    
    # Normalize flat by average of all modulations.
    if norm_method == "avg":
        for lambd in range(N_wls):
            for mod in range(N_mods):
                norm_ff[0, lambd, mod] = data[0, lambd, mod] / np.mean(data[0, lambd, :, 300:-300, 300:-300])
                norm_ff[1, lambd, mod] = data[1, lambd, mod] / np.mean(data[1, lambd, :, 300:-300, 300:-300])
    
    # Normalize flat by each modulation separately.
    elif norm_method == "mod":
        for lambd in range(N_wls):
            for mod in range(N_mods):
                norm_ff[0, lambd, mod] = data[0, lambd, mod] / np.mean(data[0, lambd, mod, 300:-300, 300:-300])
                norm_ff[1, lambd, mod] = data[1, lambd, mod] / np.mean(data[1, lambd, mod, 300:-300, 300:-300])

    else:
        raise Exception("Invalid normalization method. Please select 'avg' or 'mod'.")

    if remove_prefilter:
        if verbose:
            print(f"\nRemove prefilter activated.")
        if volts == "read":
            ff_header = flat_obs.get_info()
            volts = [ff_header["Images_headers"][f"wv_{lambd}"][f"Mod_0"]["hvps_read_volts"][0] for lambd in range(N_wls)]

            if verbose:
                print(f"Using read voltages: {volts}")
        else:
            volts = cf.om_config[om]["V_array"]
            print("Using fixed voltages.")

        if pref_model is None:
            raise Exception("Please provide a prefilter model to remove from flat-fields.")

        else:
            flats_pref_removed = pr.remove_line_from_flat_fields(norm_ff, om = om, pref_model = pref_model, volts=volts, verbose = verbose)

            if compute_blueshift_flag:
                if verbose:
                    print(f"Computing blueshift...")
                bc1, bc2 = compute_blueshift(ff_no_norm = data[:, :, :, 300:-300, 300:-300], pref_model = pref_model,
                                             volts_axis = volts, line = line, cpos = cpos, plot_flag = plot_flag, plotname = plotname)
                print(f"Flat-fields computed in {round(time.time() - tic, 3)} s.")
                return flats_pref_removed, flat_obs.get_info(), bc1, bc2

            else:
                print(f"Flat-fields computed in {round(time.time() - tic, 3)} s.")
                return flats_pref_removed, flat_obs.get_info()
    else:
        print(f"Flat-fields computed in {round(time.time() - tic, 3)} s.")
        return norm_ff, flat_obs.get_info()
    
def cog(input_data, wave_axis, cpos = -1):
    
    l, sy, sx = input_data.shape

    Ic = input_data[cpos, :, :] # Intensity of the Continuum

    t1 = wave_axis[:, np.newaxis, np.newaxis] * (Ic[np.newaxis, :, :] - input_data) # lambda * (Ic - I)
    
    tc = (Ic[np.newaxis, :, :] - input_data) # Ic - I
    
    Itn = simps(t1, x=wave_axis, axis=0) / simps(tc, x=wave_axis,axis=0) # Integral

    return Itn


def compute_blueshift(ff_no_norm, pref_model, volts_axis, line, cpos = -1, plot_flag = False, plotname = "blueshift.png"):

    def fit_2d_surface(data):
        ny, nx = data.shape
        xg, yg = np.meshgrid(np.arange(nx), np.arange(ny)) # Generate a meshgrid for the fitting

        # Flatten everything
        x_flat = xg.ravel()     # Flatten the x-coordinates
        y_flat = yg.ravel()     # Flatten the y-coordinates
        z_flat = data.ravel()   # Flatten the data

        # --- Polynomial model (3rd degree) ---
        def poly2d(X, a00, a10, a01, a20, a11, a02, a03, a30, a21, a12):
            x, y = X  
            return (
                a00
                + a10 * x
                + a01 * y
                + a20 * x**2
                + a11 * x * y
                + a02 * y**2 
                + a03 * x**3
                + a30 * y**3
                + a21 * x**2 * y
                + a12 * x * y**2
            ) 
        
        # Fit
        popt, _ = curve_fit(poly2d, (x_flat, y_flat), z_flat) # Fit to data.

        # Evaluate model
        z_fit_flat = poly2d((x_flat, y_flat), *popt)     # Generate surface
        z_fit = z_fit_flat.reshape((ny, nx))             # Reshape to 2d

        return z_fit
    
    if isinstance(pref_model, str):
        pref_model = pr.load_model(pref_model)

    shape = np.shape(ff_no_norm)
    npix = shape[-1]
    nlambda = shape[1]
    avg_flat = np.mean(ff_no_norm, axis = 2) # Average over modulations

    fitted_surface = np.zeros((2, npix, npix))
    # Define x and y for spline. 
    x = np.arange(npix)
    y = np.arange(npix)

    # FIT surface for each camera. 
    for cam in range(2):
        smoothed_flat = gaussian_filter(np.mean(avg_flat[cam], axis = 0), sigma=31) # Smooth to remove spatial variations. 
        spline_model = RectBivariateSpline(x, y, smoothed_flat, kx=2, ky=2)  # Fit a 2D surface to remove spatial variation of flat. 
        fitted_surface[cam] = spline_model(x, y)

    # Remove spatial dependence of the flat-field and prefilter effect
    flattened_flat = np.zeros_like(avg_flat)
    for cam in range(2):
        for lambd_ind, volts in enumerate(volts_axis):
            flattened_flat[cam, lambd_ind] = (avg_flat[cam, lambd_ind] / fitted_surface[cam]) / pref_model(volts) 

    central_wave = cf.prefilters_config[line]["l0"]
    wave_axis = pr.volts_2_lambda(np.array(volts_axis), cf.prefilters_config[line]) # Convert volts to wavelength.

    # Get the central wavelength
    cog_map_cam1 = cog(flattened_flat[0], wave_axis, cpos = cpos) - central_wave
    cog_map_cam2 = cog(flattened_flat[1], wave_axis, cpos = cpos) - central_wave

    # Smoothing for better fitting.
    smoothed_cog_cam1 = gaussian_filter(cog_map_cam1, sigma=35)
    smoothed_cog_cam2 = gaussian_filter(cog_map_cam2, sigma=35)
        
    # Fitting the 2d surface to data. 
    blueshift_cam1 = fit_2d_surface(smoothed_cog_cam1)
    blueshift_cam2 = fit_2d_surface(smoothed_cog_cam2)

    if plot_flag:
        print("Plotting.")
        # Auxiliary function to make colorbars
        def colorbar(mappable):
            ax = mappable.axes
            fig = ax.figure
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            return fig.colorbar(mappable, cax=cax)
        # Auxiliary function to plot blueshift
        def plot_blueshift(ax, data):
            im = ax.imshow(data, cmap = 'Spectral', origin = 'lower')
            colorbar(im)
            contour_levels = np.linspace(data.min(), data.max(), 10)  
            contours = ax.contour(data, levels=contour_levels, colors='w', linewidths=1)
            ax.clabel(contours, inline=True, fontsize=8, fmt="%.3f")

        fig, axs = plt.subplots(2, 3, figsize=(12, 8))
        im = axs[0, 0].imshow(cog_map_cam1, cmap = 'Spectral', origin = 'lower')
        colorbar(im)
        im = axs[0, 1].imshow(smoothed_cog_cam1, cmap = 'Spectral', origin = 'lower')
        colorbar(im)
        plot_blueshift(axs[0, 2], blueshift_cam1)
        im = axs[1, 0].imshow(cog_map_cam1, cmap = 'Spectral', origin = 'lower')
        colorbar(im)
        im = axs[1, 1].imshow(smoothed_cog_cam1, cmap = 'Spectral', origin = 'lower')
        colorbar(im)
        plot_blueshift(axs[1, 2], blueshift_cam1)

        fig.suptitle(f"Blueshift computation")
        axs[0, 0].set_ylabel("Camera 1")
        axs[1, 0].set_ylabel("Camera 2")
        axs[1, 0].set_xlabel("Original COG")
        axs[1, 1].set_xlabel("Smoothed COG")
        axs[1, 2].set_xlabel("Blueshift")
        fig.tight_layout()
        
        print(plotname)
        plt.savefig(plotname, bbox_inches='tight') 
        plt.show()
                    
        return blueshift_cam1, blueshift_cam2

def correct_observation(data, ff, onelambda = False):
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
        data = data[:, np.newaxis] # To allow for only one lamdba.

    # Get shape for data
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


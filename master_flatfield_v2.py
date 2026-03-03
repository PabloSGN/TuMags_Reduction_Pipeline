# ---------------------------- DESCRIPTION --------------------------------------- #

"""
Function to compute the master flat field from a set of flat-field observations of 
a single observation mode. 

Modified by FJBM on 2025-10-28.
2nd version (2026-02-24): changed zero blueshift to coincide with minimum blueshift across FOV
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Standard library
import time
from pathlib import Path

# Third-party libraries
import numpy as np
from tqdm import tqdm
from scipy.interpolate import interp1d
from scipy.optimize import minimize

# Local application modules
import config as cf
import image_handler as ih
import prefilter_removal_fjbm as pr

# ------------------------------  CODE  ------------------------------------------ # 

def compute_master_flat_field(flat_fields_paths, dc, lambda_repeat = 4, verbose = False, 
                              norm_method = "blueshift", remove_prefilter = True,
                              pref_model = True, import_blueshift_guess=True,
                              volts = None, modify_linearity = ([1539,1540],[1.0,1.0]),
                              norm_roi = [300,-300,300,-300]):

    """
    Function to compute the flat-field observation from the images paths. 

    inputs:
        - flat_field_paths (list) : List containing all the paths to the images composing a single 
        flat-field observation. 
        - dc (np.array) : Dark current. 
        - lambda_repeat (int, default : 4) : Lambda repeat parameter of the observation.
        - norm_method (str, default : "blueshift") : Normalization method. avg, mod or blueshift (FJBM).
        - remove_prefilter (Boolean, default : False) : Option to remove prefilter from the flats profiles. 
        - pref_model : Prefilter model file rerquired if remove_prefilter = True. 
        - volts (str / None, default = None) : Set to "read" if read voltages are to be used for the pref_removal.
        If None, fixed voltages are used. 
        - import_blueshift_guess (Boolean, default : False) : Option to import previous
          blueshift guess to speed up the fitting process (FJBM).
    returns:
        - ff_norm (np.array) : Array containing the flat field (cams, Nlambda, Nmods, Nx, Ny)
        - ff_info (dictionary) : Dictionary containing all info of the flat-field.

    COMMENTS: modified by FJBM on 2025-12-15 to compute blueshift maps and apply
    the flat field correction defined in De la Cruz Rodriguez et al. 2015.    
    """
    
    tic = time.time()

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
    wvlv=np.array(cf.om_config[om]["lambda_array"]) * 1e-13  # Wavelengths in meters


    # Calculate number of repetitions
    if len(flat_fields_paths) % (2 * N_wls * N_mods * lambda_repeat) == 0:
        nreps = int(len(flat_fields_paths) / (2 * N_wls * N_mods * lambda_repeat))
    else:
        raise Exception("Observations are incomplete, please remove images from incomplete OC. This will be upgraded...")
    
    if verbose:
        print(f"Observation Mode: {om}")
        print(f"Nº of repetitions: {nreps}")
        print(f"Nº of wavelengths: {N_wls}")
        print(f"Nº of Modulations: {N_mods}")
        print(f"Wavelength sampling (pm): {wvlv*1e12}")


    # Read images and correct them from dark current.
    flat_obs = ih.nominal_flat(om, flat_fields_paths, nreps, dc, modify_linearity = modify_linearity)
    data = flat_obs.get_data()

    #Define range of pixels to be used for blueshift calculation and to compute the mean profile
    x0=170
    y0=170
    xf=1800
    yf=1800
    x0_mean=1200
    y0_mean=1200
    xf_mean=x0_mean+100
    yf_mean=y0_mean+100

    #Compute mean value over central region for blueshift calculation
    cam1_ave=np.mean(data[0,:,:,x0_mean:xf_mean,y0_mean:yf_mean],axis=(1,2,3))
    cam2_ave=np.mean(data[1,:,:,x0_mean:xf_mean,y0_mean:yf_mean],axis=(1,2,3))

    norm_ff=np.ones(data.shape) #Initialize normalized flat field array

    #Normalize flat field according to selected method
    if norm_method=="blueshift":
        #Interpolation of the averaged profiles over wavelength
        cam1_interp=interp1d(wvlv,cam1_ave,kind='quadratic',bounds_error=False,
                            fill_value='extrapolate')

        cam2_interp=interp1d(wvlv,cam2_ave,kind='quadratic',bounds_error=False,
                            fill_value='extrapolate')
        
    

        """
        Fit individual profiles to mean profile to find blueshifts and
        the scaling factor
        """
        meth='Nelder-Mead' #Optimization method for scipy minimize function

        #Import guess from cam1 results of COMM1
        if import_blueshift_guess:
            if om=='2.02' or om=='2.06':
                fname_guess='ff_om2.02_D10-4340-5619_nonorm'
            elif om=='1':  
                fname_guess='ff_om1_D10-2740-4339_nonorm'  

            BASE_DIR = Path(__file__).resolve().parent
            file_path_blueshift = BASE_DIR / f"{fname_guess}_wvl_shifts_cam1.npy"
            blueshift_guess = np.load(file_path_blueshift)
            file_path_scale = BASE_DIR / f"{fname_guess}_scales_cam1.npy"
            scale_guess=np.load(file_path_scale)

        #Fit all pixels for the selected camera
        Nx=xf-x0
        Ny=yf-y0
        wvl_shifts=np.zeros((2,Nx,Ny))
        scales=np.zeros((2,Nx,Ny))

        print('Fitting blueshift and scale for each camera...')
        for cam in range(2):
            print('Camera ',cam+1)
            if cam==0:
                cam_ave=cam1_ave
            else:
                cam_ave=cam2_ave

            #Interpolation of the averaged profiles over wavelength
            if om=='2.02': #Add a "fake" point to avoid divergence when extrapolating
                wvlv_interp=np.append(wvlv,30e-12) #Fake wavelength close to continuum
                cam_ave_interp=np.append(cam_ave,cam_ave[-1]) #Fake point equal to continuum
            else: #Do nothing for other modes
                wvlv_interp=wvlv
                cam_ave_interp=cam_ave

            cam_interp=interp1d(wvlv_interp,cam_ave_interp,kind='quadratic',bounds_error=False,
                                fill_value='extrapolate')
            
            #Compute line center and set the zero blueshift to coincide with it. 
            wvlv2=np.linspace(wvlv[0],wvlv[-1],5000)
            line_center=wvlv2[np.argmin(cam_interp(wvlv2))]

            #Fit each pixel
            for i in tqdm(range(Nx)):
                for j in range(Ny):
                    #Mean spectral profile over the four modulations excluding continuum point
                    ind_fit0=0 #Initial wavelength index to fit
                    ind_fitf=-1 #Final wavelength index to fit
                    cam_indiv=np.mean(data[cam,ind_fit0:ind_fitf, :,x0+i,y0+j], axis=(1))           

                    #Define merit function to fit wavelength shift and scaling factor
                    def merit_indiv1(params):
                        delta_wvl=params[0]
                        scale=params[1]

                        #Shifted wavelength axis
                        wvlv_shifted=wvlv[ind_fit0:ind_fitf]+delta_wvl

                        #Interpolated mean profile at shifted wavelengths
                        if cam==0:
                            cam_mean_shifted=cam1_interp(wvlv_shifted)
                        else:
                            cam_mean_shifted=cam2_interp(wvlv_shifted)

                        #Compute difference with individual profile
                        diff=scale*cam_mean_shifted - cam_indiv
                        return np.sum(diff**2)
                    
                    
                    #Fit wavelength shift and scale
                    if import_blueshift_guess:
                        guess_ij=[blueshift_guess[i,j],scale_guess[i,j]]
                    else:
                        guess_ij=[0,1]
                    
                    minim=minimize(merit_indiv1,x0=guess_ij,
                                method=meth,bounds=[(-10e-12,10e-12),(0.1,2)])

                    #Save results in arrays
                    wvl_shifts[cam,i,j]=minim.x[0]
                    scales[cam,i,j]=minim.x[1]

            #Compute offset of the wavelength sampling with respect to the line center 
            max_shift=np.max(wvl_shifts[cam,:,:])
            wvl_offset=line_center+max_shift
            if verbose:
                print('Sampling offset with respect to line center for cam %g(pm): '%(cam+1),wvl_offset*1e12)

            #Correct zero blueshift to coincide with minimum blueshift across FOV
            wvl_shifts[cam,:,:]-=max_shift


            """
            Interpolate the spectral profiles of the map to correct for the blueshift
            """
            for i in tqdm(range(Nx)):
                for j in range(Ny):
                    cam_shifted_average=0
                    for k in range(4):
                        #Function with interpolation information
                        if om=='2.02': #Add point equal to continuum
                            cam_data_interp=np.append(data[cam,:,k,x0+i,y0+j],data[cam,-1,k,x0+i,y0+j]) 
                        else:
                            cam_data_interp=data[cam,:,k,x0+i,y0+j]  
                        cam_interp=interp1d(wvlv_interp,cam_data_interp,kind='quadratic',bounds_error=False,
                                    fill_value='extrapolate')

                        #Compute shifted profile averaged over the four modulations
                        cam_shifted=cam_interp(wvlv - wvl_shifts[cam,i,j])
                        cam_shifted_average+=cam_shifted
                    cam_shifted_average/=4

                    #Apply Eq. 3 from De la Cruz et al. 2015 to each modulation of the flat   
                    for k in range(4):
                        norm_ff[cam,:,k,x0+i,y0+j]=scales[cam,i,j]*data[cam,:,k,x0+i,y0+j]/cam_shifted_average
    #Normalize flat by average of all modulations.
    elif norm_method == "avg":
        norma = np.mean(data[:, :, :, norm_roi[0]:norm_roi[1], norm_roi[2]:norm_roi[3]],axis=(2,3,4))
        for lambd in range(N_wls):
            for mod in range(N_mods):
                norm_ff[0, lambd, mod] = data[0, lambd, mod] / norma[0, lambd]
                norm_ff[1, lambd, mod] = data[1, lambd, mod] / norma[1, lambd]
    
    #Normalize flat by each modulation separately.
    elif norm_method == "mod":       
        norma = np.zeros(np.shape(data[:,:,:]))
        for lambd in range(N_wls):
            for mod in range(N_mods):
                norma[:,lambd,mod] = np.mean(data[:, lambd, mod, norm_roi[0]:norm_roi[1], norm_roi[2]:norm_roi[3]],axis=(1,2))
                norm_ff[0, lambd, mod] = data[0, lambd, mod] / norma[0,lambd,mod]
                norm_ff[1, lambd, mod] = data[1, lambd, mod] / norma[1,lambd,mod]
    elif norm_method == "none" or norm_method == None:
        norm_ff = data
    else:
        raise Exception("Invalid normalization method. Please select 'blueshift','avg', 'mod' or 'none'")


    """
    Remove prefilter if required
    """
    if remove_prefilter:
        if pref_model is True:
            prefilter=pr.prefilter_model(om,wvlv) 
        else:
            prefilter=pr.prefilter_fitting(cam1_ave,om,wvlv) #Use cam1 for fitting
        norm_ff*=prefilter[np.newaxis,:,np.newaxis,np.newaxis,np.newaxis]
        
    print(f"Flat-fields computed in {round(time.time() - tic, 3)} s.")
    return norm_ff, flat_obs.get_info()
    
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


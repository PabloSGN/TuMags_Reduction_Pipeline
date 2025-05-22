# ---------------------------- DESCRIPTION --------------------------------------- #

"""
Function to compute the master flat field from a set of flat-field observations of 
a single observation mode. 
"""

# ------------------------------ IMPORTS ----------------------------------------- #
import numpy as np
import time

# Own modules
from . import config as cf
from . import image_handler as ih
from . import prefilter_removal as pr

# ------------------------------  CODE  ------------------------------------------ # 

def compute_master_flat_field(flat_fields_paths, dc, lambda_repeat = 4, verbose = False, 
                              norm_method = "avg", remove_prefilter = False, pref_model = None, 
                              volts = None):
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
        if volts == "read":
            ff_header = flat_obs.get_info()
            volts = [ff_header["Images_headers"][f"wv_{lambd}"][f"Mod_0"]["hvps_read_volts"][0] for lambd in range(N_wls)]
            print(volts)

            if verbose:
                print(f"Using read voltages: {volts}")
        else:
            print("Using fixed voltages.")

        if pref_model is None:
            raise Exception("Please provide a prefilter model to remove from flat-fields.")
        else:
            flats_pref_removed = pr.remove_line_from_flat_fields(norm_ff, om = om, pref_model = pref_model, volts=volts, verbose = verbose)

            print(f"Flat-fields computed in {round(time.time() - tic, 3)} s.")
            return flats_pref_removed, flat_obs.get_info()
    else:
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


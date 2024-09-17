# ---------------------------- DESCRIPTION --------------------------------------- #

"""
Function to compute the master flat field from a set of flat-field observations of 
a single observation mode. 
"""

# ------------------------------ IMPORTS ----------------------------------------- #
import numpy as np
import time

import config as cf
import image_handler as ih

# ------------------------------ CONFIG ------------------------------------------ #

# ------------------------------  CODE  ------------------------------------------ # 

def compute_master_flat_field(flat_fields_paths, dc, lambda_repeat = 4, verbose = False):
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
    naccs = h["nAcc"]

    if verbose:
        print(f"Observation Mode: {om}")
        print(f"Nº of repetitions: {nreps}")
        print(f"Nº of wavelengths: {N_wls}")
        print(f"Nº of Modulations: {N_mods}")

    # Read images and correct them from dark current (scaled to accumulation number)
    flat_obs = ih.nominal_flat(om, flat_fields_paths, nreps, dc * naccs)

    data = flat_obs.get_data()

    norm_ff = np.zeros(np.shape(data))
    # Normalize flat-fields
    for lambd in range(N_wls):
        for mod in range(N_mods):
            norm_ff[0, lambd, mod] = data[0, lambd, mod] / np.mean(data[0, lambd, mod, 300:-300, 300:-300])
            norm_ff[1, lambd, mod] = data[1, lambd, mod] / np.mean(data[1, lambd, mod, 300:-300, 300:-300])

    print(f"Flat-fields computed in {round(time.time() - tic, 3)} s.")
    return norm_ff, flat_obs.get_info()
    




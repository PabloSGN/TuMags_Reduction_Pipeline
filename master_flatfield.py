# ---------------------------- DESCRIPTION --------------------------------------- #

# ------------------------------ IMPORTS ----------------------------------------- #
import numpy as np
import time

from utils import read_Tumag
import config as cf
import image_handler as ih

# ------------------------------ CONFIG ------------------------------------------ #

# ------------------------------  CODE  ------------------------------------------ # 

def compute_master_flat_field(flat_fields_paths, dc, lambda_repeat = 4, verbose = False):
    tic = time.time()

    if verbose:
        print(f"Computing flats.")
        print(f"N flats: {len(flat_fields_paths)}")

    # Read first image to get info from header
    _, h = ih.read(flat_fields_paths[0])

    om = h["ObservationMode"]
    N_wls = cf.om_config[om]["Nlambda"]
    N_mods = cf.om_config[om]["Nmods"]
    nreps = int(len(flat_fields_paths) / (2 * N_wls * N_mods * lambda_repeat))

    if verbose:
        print(f"Observation Mode: {om}")
        print(f"Nº of repetitions: {nreps}")
        print(f"Nº of wavelengths: {N_wls}")
        print(f"Nº of Modulations: {N_mods}")

    # Read images
    flat_obs = ih.nominal_flat(om, flat_fields_paths, nreps, dc)

    print(f"Flat-fields computed in {round(time.time() - tic, 3)} s.")
    return flat_obs.get_data(), flat_obs.get_info()
    




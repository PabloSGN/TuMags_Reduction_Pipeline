# ---------------------------- DESCRIPTION --------------------------------------- #
"""
Function to compute the dark current from a set of observations. 

It averages all observations and retruns the dark current per acumulation.

"""
# ------------------------------ IMPORTS ----------------------------------------- #

import numpy as np
from utils import read_Tumag
import time

# ------------------------------  CODE  ------------------------------------------ # 

def compute_master_darks(dark_paths, verbose = False):
    """
    Function to compute the dark current from the paths to the images. 
    inputs:
        - dark_paths (list) : List of paths to the dark current images.
    returns: 
        - dc (np.array) : Array of the dark current for both cameras.
    
    """

    tic = time.time()
    # Read first image to obtain image size.
    first_dark, head = read_Tumag(dark_paths[0]) 
    dark_current = np.zeros((2, np.shape(first_dark)[0],np.shape(first_dark)[0]))
    
    darks_cam_1 = [x for x in dark_paths if "_0_" in x]
    darks_cam_2 = [x for x in dark_paths if "_1_" in x]

    # Proccessing cam 1 darks
    if verbose:
        print(f"\nComputing darks.")
        print(f"------------------")
        print(f"N darks for cam1 : {len(darks_cam_1)}")
        print(f"N darks for cam2 : {len(darks_cam_2)}")
        print(f"N accumulations : {head['nAcc']}")


    for _, img_path in enumerate(darks_cam_1):
        I, _ = read_Tumag(img_path)
        dark_current[0] += I

    for _, img_path in enumerate(darks_cam_2):
        I, _ = read_Tumag(img_path)
        dark_current[1] += np.flip(I, axis = -1)

    dark_current[0] /= len(darks_cam_1)
    dark_current[1] /= len(darks_cam_2)

    # We compute the dark current per accumulation to scale it to other observations.
    dark_current /= head['nAcc']

    print(f"Dark current computed in {round(time.time() - tic, 3)} s.\n")
    return dark_current







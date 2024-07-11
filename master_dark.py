# ---------------------------- DESCRIPTION --------------------------------------- #

# ------------------------------ IMPORTS ----------------------------------------- #
import glob
import os
import numpy as np
from utils import read_Tumag

import time
import matplotlib.pyplot as plt

from field_stop_finder import find_fieldstop

# ------------------------------ CONFIG ------------------------------------------ #

# ------------------------------  CODE  ------------------------------------------ # 

def compute_master_darks(darks_cam_1, darks_cam_2, fieldstop_cam1, fieldstop_cam2, verbose = False):

    tic = time.time()
    # Read first image to obtain image size.
    first_dark, _ = read_Tumag(darks_cam_1[0]) 
    dark_current = np.zeros((2, np.shape(first_dark)[0],np.shape(first_dark)[0]))
    
    # Proccessing cam 1 darks
    if verbose:
        print(f"Computing darks ...")
        print(f"N darks for cam1 : {len(darks_cam_1)}")
        print(f"N darks for cam2 : {len(darks_cam_1)}")

    for img_ind, img_path in enumerate(darks_cam_1):
        I, _ = read_Tumag(img_path)
        dark_current[0] += I

    for img_ind, img_path in enumerate(darks_cam_2):
        I, _ = read_Tumag(img_path)
        dark_current[1] += I
        
    field_stop_mask_cam1 = np.zeros((np.shape(first_dark)))
    field_stop_mask_cam2 = np.zeros((np.shape(first_dark)))

    field_stop_mask_cam1[fieldstop_cam1[0][0] : fieldstop_cam1[0][1], fieldstop_cam1[1][0] : fieldstop_cam1[1][1]] = 1
    field_stop_mask_cam2[fieldstop_cam2[0][0] : fieldstop_cam2[0][1], fieldstop_cam2[1][0] : fieldstop_cam2[1][1]] = 1

    dark_current[0] /= len(darks_cam_1)
    dark_current[1] /= len(darks_cam_2)

    dark_current[0] *= field_stop_mask_cam1
    dark_current[1] *= field_stop_mask_cam2

    print(f"Dark current computed in {round(time.time() - tic, 3)} s.")
    return dark_current







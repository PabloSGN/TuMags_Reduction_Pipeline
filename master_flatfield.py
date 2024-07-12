# ---------------------------- DESCRIPTION --------------------------------------- #

# ------------------------------ IMPORTS ----------------------------------------- #
import numpy as np
import time

from utils import read_Tumag

# ------------------------------ CONFIG ------------------------------------------ #

# ------------------------------  CODE  ------------------------------------------ # 

def compute_master_flat_field(flat_fields_paths_cam1, flat_fields_paths_cam2,
                              dc, fieldstop_cam1, fieldstop_cam2,
                              verbose = False):
    tic = time.time()
    
    # standard flat field measurements 
    first_flat, _ = read_Tumag(flat_fields_paths_cam1[0]) 
    ff = np.zeros((2, np.shape(first_flat)[0],np.shape(first_flat)[0]))

    if verbose:
        print(f"Computing flats ...")
        print(f"N flats for cam1 : {len(flat_fields_paths_cam1)}")
        print(f"N flats for cam2 : {len(flat_fields_paths_cam2)}")

    for img_ind, img_path in enumerate(flat_fields_paths_cam1):
        I, _ = read_Tumag(img_path)
        ff[0] += I - dc[0]

    for img_ind, img_path in enumerate(flat_fields_paths_cam2):
        I, _ = read_Tumag(img_path)
        I = np.flip(I, axis = -1)
        ff[1] += I - dc[1]
        
    ff[0] /= len(ff)
    ff[1] /= len(ff)

    field_stop_mask_cam1 = np.zeros((np.shape(first_flat)))
    field_stop_mask_cam2 = np.zeros((np.shape(first_flat)))
    field_stop_mask_cam1[fieldstop_cam1[0][0] : fieldstop_cam1[0][1], fieldstop_cam1[1][0] : fieldstop_cam1[1][1]] = 1
    field_stop_mask_cam2[fieldstop_cam2[0][0] : fieldstop_cam2[0][1], fieldstop_cam2[1][0] : fieldstop_cam2[1][1]] = 1

    ff[0] *= field_stop_mask_cam1
    ff[1] *= field_stop_mask_cam2

    print(f"Flat-fields computed in {round(time.time() - tic, 3)} s.")
    return ff
    




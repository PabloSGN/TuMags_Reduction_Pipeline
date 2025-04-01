# ---------------------------- DESCRIPTION --------------------------------------- #

"""
Phase diversity module. Parses the images, and computes PD from raw data. 

authors: 
Pablo Santamarina Guerrero(pablosantamarinag@gmail.com)
Francisco Javier Bailén (fbailen@iaa.es)

Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in Libs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time

# Own libs
import config as cf
import image_handler as ih
from master_dark import compute_master_darks
from master_flatfield import compute_master_flat_field
import pd_functions_v22 as pdf

# ------------------------------  CODE  ------------------------------------------ # 

def pd_observation_parser(indexes, verbose = False):

    """
    Function that given the indexes of a PD observation reorders the paths and
    returns a numpy array with the ordered paths. 

    Inputs: 
        - indexes (str) : query for the image_handler.get_images_paths() function.
                    format must be: DXX-YYYY
        - verbose (bool) : Boolean to print info of the images. Default : False
    Returns:
        - images_ordered (np.array) : Numpy array with the images paths ordered.
                                      Order [camera - Filter - PD/F4 - Nimages]
        - Nims_pd (int) : Number of images per camera filter and PD/F4
    
    """
    image_paths = ih.get_images_paths(indexes)

    if len(image_paths) == 320:
        Nims_pd = 40
    elif len(image_paths) == 256:
        Nims_pd = 32
    else:
        raise Exception(f"Error in PD processing. Expected 320 or 256 images but got {len()}.")

    try:
        # Reorder array: [Filter - PD/F4 - Nimages - Cameras]
        images_ordered = np.array(image_paths).reshape(2, 2, Nims_pd, 2) 
    except:
        raise Exception(f"\nIncorrect number of images for pd modes.\n{len(image_paths)} Images found but 320 were expected (2 (filters) x 2 (PD and F4) x 2 (cams) x 4 (Nimages per mode)).")

    # Change array order: [cam, filter, PD/F4, Nimage]
    images_ordered = np.transpose(images_ordered, axes = (3, 0, 1, 2)) 

    if verbose:
        print(f"\nParsing PD images....")
        print(f"------------------")
        print(f"Total number of images : {len(image_paths)}")
        
        _, H1 =ih.read(images_ordered[0, 0, 0, 0])
        _, H2 =ih.read(images_ordered[0, 1, 0, 0])
        
        print(f"Filter 0 : {H1['FW2']}")
        print(f"Filter 1 : {H2['FW2']}")

    return images_ordered, Nims_pd

def process_pd_observation(indexes, filt, verbose = False):

    print("Processing PD observation...")

    paths, Nims_pd = pd_observation_parser(indexes, verbose = verbose)

    paths = paths[:, filt] # Select filter to process. 

    images = np.zeros((2, 2, Nims_pd, cf.xsize, cf.ysize))
    # Loop over pd and f4 images
    for pd in range(2):
        # Loop over all images.
        for Nimg in range(Nims_pd):
            I, _ = ih.read(paths[0, pd, Nimg])
            images[0, pd, Nimg] = I

            # Flip cam 2 to match flats and dc
            I, _ = ih.read(paths[1, pd, Nimg])
            images[1, pd, Nimg] = np.flip(I, axis = -1)

    return images

def apply_reconstruction(data, ID = None, zkes = None, verbose = True):
    """
    Function that applies the restore ima routine to all images of an observing
    mode. 

    Inputs: 
        - data : Numpy array containing obs mode.
        - ID : Sunrise ID for automatic zernike identification 
        - zkes : If direct zernikes are passed.
    Returns: 
        - Reconstructed data
    """
    tic = time.time()

    if zkes is None:
        zkes = pdf.import_zernikes(ID)

    # Get shape for data
    shape = np.shape(data)
    nlambda = shape[1]
    nmods = shape[2]

    reconstructed = np.zeros_like(data)

    print("Applying pd reconstruction...")
    
    for lamb in range(nlambda):
        for mod in range(nmods):
            if verbose:
                print(f"Wl - {lamb} - Mod - {mod}")
            for cam in range(2):
                reconstructed[cam, lamb, mod], _ = pdf.restore_ima(data[cam, lamb, mod], zkes)

    tac = time.time()
    if verbose:
        print(f"Wavefront reconstruction finished in {round(tac - tic, 3)}s.")
    print("Reconstruction finished...")

    return reconstructed, zkes




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

# Own libs
import config as cf
import image_handler as ih
from master_dark import compute_master_darks
from master_flatfield import compute_master_flat_field

# Config

# At the start of the mission
Nims_pd = 40 # Images of PD per camera

# from the middle of the mission
#Nims_pd = 32 # Images of PD per camera
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
    
    """
    image_paths = ih.get_images_paths(indexes)

    try:
        # Reorder array: [Filter - PD/F4 - Nimages - Cameras]
        images_ordered = np.array(image_paths).reshape(2, 2, 40, 2) 
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

    return images_ordered

def process_pd_observation(indexes, filt, verbose = False):

    print("Processing PD observation...")

    paths = pd_observation_parser(indexes, verbose = verbose)

    paths = paths[:, filt] # Select filter to process. 

    images = np.zeros((2, 2, Nims_pd, cf.xsize, cf.ysize))
    # Loop over 2 cameras
    for cam in range(2):
        # Loop over pd and f4 images
        for pd in range(2):
            # Loop over all images.
            for Nimg in range(Nims_pd):
                print(Nimg)
                I, _ = ih.read(paths[cam, pd, Nimg])
                images[cam, pd, Nimg] = I

    return images
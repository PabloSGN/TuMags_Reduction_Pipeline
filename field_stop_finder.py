# ---------------------------- DESCRIPTION --------------------------------------- #

"""
author: Pablo Santamarina Guerrero(pablosantamarinag@gmail.com) 
Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in Libs
import glob
import numpy as np
import time

from scipy.signal import correlate, correlate2d

import matplotlib.pyplot as plt

# Own libs
from utils import read_Tumag
import config as cf

# ------------------------------- CONFIG ----------------------------------------- #

# Central pinhole
central_pinhole_x =  [875, 1025]
central_pinhole_y =  [925, 1075]

# Small pinhole
# central_pinhole_x = [400, 500]
# central_pinhole_y = [1300, 1400]

margin = 10 # Margin from fieldstop

# ------------------------------  CODE  ------------------------------------------ # 

def find_fieldstop(cam1 = None, verbose = False, plot_flag = False, margin = margin):

    tic = time.time()

    if verbose:
        print("Finding fieldstop field stop...")
    
    if plot_flag:
        fig, axs  = plt.subplots(1, 2,figsize = (10, 5))
        axs[0].imshow(cam1, origin = 'lower', cmap = 'gray')
    
    size = np.shape(cam1)[0]

    # Position to find cuts
    lines = np.linspace(0, size, 7)
    lines = [int(x) for x in lines[1:-1]]
    
    # Looking for cuts
    hcuts_left_c1 = []
    hcuts_right_c1 = []
    vcuts_top_c1 = []
    vcuts_bottom_c1 = []
    for l in lines:
        # Camera 1
        hcut1 = np.argmax(np.gradient(cam1[l, :]))
        hcut2 = np.argmin(np.gradient(cam1[l, :]))
        vcut1 = np.argmax(np.gradient(cam1[:, l]))
        vcut2 = np.argmin(np.gradient(cam1[:, l]))
        hcuts_right_c1.append(hcut1)
        hcuts_left_c1.append(hcut2)                                    
        vcuts_top_c1.append(vcut1)
        vcuts_bottom_c1.append(vcut2)

        if plot_flag:
            axs[0].plot([l, l], [0, size], color = 'crimson' , lw = 1)
            axs[0].plot([0, size], [l, l], color = 'crimson' , lw = 1)
            axs[0].scatter(l, vcut1, marker = 'x', c = 'dodgerblue')
            axs[0].scatter(l, vcut2, marker = 'x', c = 'darkorange')
            axs[0].scatter(hcut1, l, marker = 'x', c = 'dodgerblue')
            axs[0].scatter(hcut2, l, marker = 'x', c = 'darkorange')

    # Selecting the innermost points (in case border is tilted)

    print(vcuts_bottom_c1)
    vcut_right_c1 = np.min(hcuts_left_c1) - margin 
    vcut_left_c1 = np.max(hcuts_right_c1) + margin
    hcut_top_c1 = np.min(vcuts_bottom_c1) - margin
    hcut_bottom_c1 = np.max(vcuts_top_c1) + margin

    cam1_fieldstop = np.array([[hcut_bottom_c1, hcut_top_c1], [vcut_left_c1, vcut_right_c1]])

    if plot_flag:
        axs[0].plot([vcut_right_c1, vcut_right_c1], [0, size], c = 'deeppink')
        axs[0].plot([vcut_left_c1, vcut_left_c1], [0, size], c = 'deeppink')
        axs[0].plot([0, size], [hcut_top_c1, hcut_top_c1], c = 'deeppink')
        axs[0].plot([0, size], [hcut_bottom_c1, hcut_bottom_c1], c = 'deeppink')
        axs[1].imshow(cam1[hcut_bottom_c1:hcut_top_c1, vcut_left_c1:vcut_right_c1], origin = 'lower', cmap = 'gray')
        axs[0].set_xlim(0, size)
        axs[0].set_ylim(0, size)
        axs[0].set_ylabel("Cam 1")
        plt.tight_layout()
        plt.show()

    if verbose:
        print(f"Cam 1 Field stop: {cam1_fieldstop}\nHeight : {cam1_fieldstop[0][1] - cam1_fieldstop[0][0]}\nWidth : {cam1_fieldstop[1][1] - cam1_fieldstop[1][0]}")

    print(f"Field stop computation finished in {round(time.time() - tic, 3)}s.")

    return cam1_fieldstop

def apply_fieldstop_and_align_single_image(cam1, cam2, field_stop_c1, field_stop_c2):

    # New array to store fieldstopped and align images
    new_cam1 = np.zeros((np.shape(cam1)))
    new_cam2 = np.zeros((np.shape(cam1)))

    # Field stop of cam 1
    field_stop_mask_cam1 = np.zeros((np.shape(cam1)))
    field_stop_mask_cam1[field_stop_c1[0][0] : field_stop_c1[0][1], field_stop_c1[1][0] : field_stop_c1[1][1]] = 1

    # Apply cam 1 field stop
    new_cam1 = field_stop_mask_cam1 * cam1

    # Align cam 2
    new_cam2[field_stop_c1[0][0] : field_stop_c1[0][1], field_stop_c1[1][0] : field_stop_c1[1][1]] = cam2[field_stop_c2[0][0] : field_stop_c2[0][1], field_stop_c2[1][0] : field_stop_c2[1][1]]    

    return new_cam1, new_cam2

def apply_fieldstop_and_align_array(data, fs_c1, fs_c2):

    shp = np.shape(data)

    # Field stop of cam 1
    field_stop_mask_cam1 = np.zeros((shp[-2:]))
    field_stop_mask_cam1[fs_c1[0][0] : fs_c1[0][1], fs_c1[1][0] : fs_c1[1][1]] = 1

    # New array to store fieldstopped and align images
    new_data = np.zeros((np.shape(data)))

    for lambd in range(shp[1]):
        for mod in range(shp[2]):
            # Apply cam 1 field stop
            new_data[0, lambd, mod] = field_stop_mask_cam1 * data[0, lambd, mod]

            # Align cam 2
            new_data[1, lambd, mod, fs_c1[0][0] : fs_c1[0][1],
                                    fs_c1[1][0] : fs_c1[1][1]] = \
                                    data[1, lambd, mod, fs_c2[0][0] : fs_c2[0][1], 
                                                        fs_c2[1][0] : fs_c2[1][1]]    

    return new_data

def compute_alignment(flat_cam1, flat_cam2 = None, pinhole_c1_path = None,
                       pinhole_c2_path = None, method = "pinhole", verbose = True, plot_flag = True, cent_ph_x = central_pinhole_x, cent_ph_y = central_pinhole_y):
    
    if method == "pinhole":
        # Read Pinhole images

        if verbose:
            print("Computing alignment with pinholes..")

        if pinhole_c1_path == None or pinhole_c2_path == None:
            raise Exception("Provide pinhole paths for pinhole alignment.")

        ph1, _ = read_Tumag(pinhole_c1_path)
        ph2, _ = read_Tumag(pinhole_c2_path)
        ph2 = np.flip(ph2, axis = -1) # Flip cam 2

        # Normalize pinholes
        ph1  = ph1 / np.max(ph1)
        ph2  = ph2 / np.max(ph2)

        # Correlation of central pinhole between 2 cams 
        if verbose:
            print("Computing correlation...")
        
        # Computing correlation of the central pinhole (normalized) between cameras. 
        correlation = correlate(ph1[cent_ph_x[0]:cent_ph_x[1], cent_ph_y[0]:cent_ph_y[1]] / np.max(ph1[cent_ph_x[0]:cent_ph_x[1], cent_ph_y[0]:cent_ph_y[1]]), 
                                ph2[cent_ph_x[0]:cent_ph_x[1], cent_ph_y[0]:cent_ph_y[1]] / np.max(ph2[cent_ph_x[0]:cent_ph_x[1], cent_ph_y[0]:cent_ph_y[1]]), mode='full', method = 'fft')
        
        # Computing shift
        max_index = np.unravel_index(np.argmax(correlation), correlation.shape) # Highest corr value.

        shift_x = max_index[0] - (cent_ph_x[1] - cent_ph_x[0] - 1) # -1 is for even sizes
        shift_y = max_index[1] - (cent_ph_y[1] - cent_ph_y[0] - 1) # -1 is for even sizes

        # Normalize flat_field and calculate field-stop
        flat_cam1 = flat_cam1 / np.max(flat_cam1)
        fs_c1 = find_fieldstop(cam1 = flat_cam1, plot_flag = plot_flag, verbose = verbose)

        # Generate field stop for cam 2 from field stop of cam 1 and shift
        fs_c2 = np.zeros(np.shape(fs_c1)).astype("int")
        fs_c2[0] = fs_c1[0] - int(shift_x)
        fs_c2[1] = fs_c1[1] - int(shift_y)

        if plot_flag:
            
            _, axs = plt.subplots(2, 2, figsize = (10, 10))
            axs[0, 0].set_title("Pinholes Diff. pre-aligned")
            im = axs[0, 0].imshow(ph1[cent_ph_x[0]:cent_ph_x[1], cent_ph_y[0]:cent_ph_y[1]] - \
                                  ph2[cent_ph_x[0]:cent_ph_x[1], cent_ph_y[0]:cent_ph_y[1]], cmap = 'afmhot')
            plt.colorbar(im, fraction=0.046, pad=0.04)

            axs[1, 0].set_title("Falts Diff. pre-aligned")
            im = axs[1, 0].imshow(flat_cam1 - flat_cam2, cmap = 'afmhot')
            plt.colorbar(im, fraction=0.046, pad=0.04)

            ph1, ph2 = apply_fieldstop_and_align_single_image(ph1, ph2, fs_c1, fs_c2)
            flat_cam1, flat_cam2_shifted = apply_fieldstop_and_align_single_image(flat_cam1, flat_cam2, fs_c1, fs_c2)

            axs[0, 1].set_title("Pinholes Diff. post-aligned")
            im = axs[0, 1].imshow(ph1[cent_ph_x[0]:cent_ph_x[1], cent_ph_y[0]:cent_ph_y[1]] - \
                                  ph2[cent_ph_x[0]:cent_ph_x[1], cent_ph_y[0]:cent_ph_y[1]], cmap = 'afmhot')
            plt.colorbar(im, fraction=0.046, pad=0.04)

            axs[1, 1].set_title("Falts Diff. post-aligned")
            im = axs[1, 1].imshow(flat_cam1 - flat_cam2_shifted, cmap = 'afmhot')
            plt.colorbar(im, fraction=0.046, pad=0.04)

            plt.tight_layout()
            plt.show()

        if verbose:
            print(f"Shift between cameras : X =  {shift_x} - Y - {shift_y} ")

        return fs_c1, fs_c2
    
    elif method == 'flats':

        if verbose:
            print("Computing alignment with flat-fields...")

        if verbose:
            print("Computing correlation...")
        # Correlation with flat-fields
        correlation = correlate(flat_cam1, flat_cam2, mode='full', method = 'fft')
        
        # Computing shift
        max_index = np.unravel_index(np.argmax(correlation), correlation.shape) # Highest corr value.

        shift_x = max_index[0] - (2016 - 1) # -1 is for even sizes
        shift_y = max_index[1] - (2016 - 1) # -1 is for even sizes

        # Normalize flat_field and calculate field-stop
        flat_cam1 = flat_cam1 / np.max(flat_cam1)
        fs_c1 = find_fieldstop(cam1 = flat_cam1, plot_flag = plot_flag, verbose = verbose)

        # Generate field stop for cam 2 from field stop of cam 1 and shift
        fs_c2 = np.zeros(np.shape(fs_c1)).astype("int")
        fs_c2[0] = fs_c1[0] - int(shift_x)
        fs_c2[1] = fs_c1[1] - int(shift_y)

        if plot_flag:
            
            _, axs = plt.subplots(1, 3, figsize = (15, 5))

            axs[0].set_title("Falts Diff. pre-aligned")
            im = axs[0].imshow(flat_cam1 - flat_cam2, cmap = 'afmhot')
            plt.colorbar(im, fraction=0.046, pad=0.04)

            flat_cam1, flat_cam2_shifted = apply_fieldstop_and_align_single_image(flat_cam1, flat_cam2, fs_c1, fs_c2)

            axs[1].set_title("Falts Diff between shifts")
            im = axs[1].imshow(flat_cam2 - flat_cam2_shifted, cmap = 'afmhot')
            plt.colorbar(im, fraction=0.046, pad=0.04)


            axs[2].set_title("Falts Diff. post-aligned")
            im = axs[2].imshow(flat_cam1 - flat_cam2_shifted, cmap = 'afmhot')
            plt.colorbar(im, fraction=0.046, pad=0.04)

            plt.tight_layout()
            plt.show()

        if verbose:
            print(f"Shift between cameras : X =  {shift_x} - Y - {shift_y} ")

        return fs_c1, fs_c2

    else:
        raise Exception("Please provide method of alignment computation : 'pinhole'")
    















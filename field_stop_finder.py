# ---------------------------- DESCRIPTION --------------------------------------- #

"""
author: Pablo Santamarina Guerrero(pablosantamarinag@gmail.com) 
Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in Libs
from astropy.io import fits
import numpy as np
import time

import matplotlib.pyplot as plt

# Own libs
from utils import read_Tumag
import config as cf

# ------------------------------  CODE  ------------------------------------------ # 

def find_fieldstop(mode = "fixed", cam1 = None, cam2 = None, verbose = False, plot_flag = False):

    if mode == "fixed":
        print("Fixed Field Stop selected.")

        c1_fieldstop = [[135, 1800], [161, 1826]]
        c2_fieldstop = [[135, 1800], [189, 1854]]

        size = np.shape(cam1)[0]

        hcut_bottom_c1 = 135
        hcut_bottom_c2 = 135
        hcut_top_c1 = 1800
        hcut_top_c2 = 1800
        vcut_left_c1 = 161 
        vcut_left_c2 = 189
        vcut_right_c1 = 1826
        vcut_right_c2 = 1854

        if plot_flag:
            fig, axs  = plt.subplots(2, 2,figsize = (10, 10))
            axs[0, 0].imshow(cam1, origin = 'lower', cmap = 'gray')
            axs[1, 0].imshow(cam2, origin = 'lower', cmap = 'gray')
            axs[0, 0].plot([vcut_right_c1, vcut_right_c1], [0, size], c = 'deeppink')
            axs[0, 0].plot([vcut_left_c1, vcut_left_c1], [0, size], c = 'deeppink')
            axs[0, 0].plot([0, size], [hcut_top_c1, hcut_top_c1], c = 'deeppink')
            axs[0, 0].plot([0, size], [hcut_bottom_c1, hcut_bottom_c1], c = 'deeppink')
            axs[1, 0].plot([vcut_right_c2, vcut_right_c2], [0, size], c = 'deeppink')
            axs[1, 0].plot([vcut_left_c2, vcut_left_c2], [0, size], c = 'deeppink')
            axs[1, 0].plot([0, size], [hcut_top_c2, hcut_top_c2], c = 'deeppink')
            axs[1, 0].plot([0, size], [hcut_bottom_c2, hcut_bottom_c2], c = 'deeppink')
            axs[0, 1].imshow(cam1[hcut_bottom_c1:hcut_top_c1, vcut_left_c1:vcut_right_c1], origin = 'lower', cmap = 'gray')
            axs[1, 1].imshow(cam2[hcut_bottom_c2:hcut_top_c2, vcut_left_c2:vcut_right_c2], origin = 'lower', cmap = 'gray')
            axs[0, 0].set_xlim(0, size)
            axs[0, 0].set_ylim(0, size)
            axs[1, 0].set_xlim(0, size)
            axs[1, 0].set_ylim(0, size)    
            axs[0, 0].set_ylabel("Cam 1")
            axs[1, 0].set_ylabel("Cam 2")

            plt.tight_layout()
            plt.show()
        if verbose:
            print(f"Cam 1 Field stop: {c1_fieldstop}")
            print(f"Cam 2 Field stop: {c2_fieldstop}")

        return c1_fieldstop, c2_fieldstop 
    
    elif mode == "auto":

        tic = time.time()

        print("Automatic field stop selected...")

        size = np.shape(cam1)[0]

        # Position to find cuts
        lines = np.linspace(0, size, 7)
        lines = [int(x) for x in lines[1:-1]]

        # Plotting image
        if plot_flag:
            fig, axs  = plt.subplots(2, 2,figsize = (10, 10))
            axs[0, 0].imshow(cam1, origin = 'lower', cmap = 'gray')
            axs[1, 0].imshow(cam2, origin = 'lower', cmap = 'gray')

        # Looking for cuts
        hcuts_left_c1 = []
        hcuts_right_c1 = []
        vcuts_top_c1 = []
        vcuts_bottom_c1 = []
        hcuts_left_c2 = []
        hcuts_right_c2 = []
        vcuts_top_c2 = []
        vcuts_bottom_c2 = []
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
                axs[0, 0].plot([l, l], [0, size], color = 'crimson' , lw = 1)
                axs[0, 0].plot([0, size], [l, l], color = 'crimson' , lw = 1)
                axs[0, 0].scatter(l, vcut1, marker = 'x', c = 'dodgerblue')
                axs[0, 0].scatter(l, vcut2, marker = 'x', c = 'darkorange')
                axs[0, 0].scatter(hcut1, l, marker = 'x', c = 'dodgerblue')
                axs[0, 0].scatter(hcut2, l, marker = 'x', c = 'darkorange')

            # Camera 2
            hcut1 = np.argmax(np.gradient(cam2[l, :]))
            hcut2 = np.argmin(np.gradient(cam2[l, :]))
            vcut1 = np.argmax(np.gradient(cam2[:, l]))
            vcut2 = np.argmin(np.gradient(cam2[:, l]))
            hcuts_right_c2.append(hcut1)
            hcuts_left_c2.append(hcut2)                                    
            vcuts_top_c2.append(vcut1)
            vcuts_bottom_c2.append(vcut2)

            if plot_flag:
                axs[1, 0].plot([l, l], [0, size], color = 'crimson' , lw = 1)
                axs[1, 0].plot([0, size], [l, l], color = 'crimson' , lw = 1)
                axs[1, 0].scatter(l, vcut1, marker = 'x', c = 'dodgerblue')
                axs[1, 0].scatter(l, vcut2, marker = 'x', c = 'darkorange')
                axs[1, 0].scatter(hcut1, l, marker = 'x', c = 'dodgerblue')
                axs[1, 0].scatter(hcut2, l, marker = 'x', c = 'darkorange')

        # Selecting the innermost points (in case border is tilted)
        vcut_right_c1 = np.min(hcuts_left_c1)
        vcut_left_c1 = np.max(hcuts_right_c1)
        hcut_top_c1 = np.min(vcuts_bottom_c1)
        hcut_bottom_c1 = np.max(vcuts_top_c1)
        vcut_right_c2 = np.min(hcuts_left_c2)
        vcut_left_c2 = np.max(hcuts_right_c2)
        hcut_top_c2 = np.min(vcuts_bottom_c2)
        hcut_bottom_c2 = np.max(vcuts_top_c2)

        cam1_fieldstop = [[hcut_bottom_c1, hcut_top_c1], [vcut_left_c1, vcut_right_c1]]
        cam2_fieldstop = [[hcut_bottom_c2, hcut_top_c2], [vcut_left_c2, vcut_right_c2]]

        if plot_flag:
            axs[0, 0].plot([vcut_right_c1, vcut_right_c1], [0, size], c = 'deeppink')
            axs[0, 0].plot([vcut_left_c1, vcut_left_c1], [0, size], c = 'deeppink')
            axs[0, 0].plot([0, size], [hcut_top_c1, hcut_top_c1], c = 'deeppink')
            axs[0, 0].plot([0, size], [hcut_bottom_c1, hcut_bottom_c1], c = 'deeppink')
            axs[1, 0].plot([vcut_right_c2, vcut_right_c2], [0, size], c = 'deeppink')
            axs[1, 0].plot([vcut_left_c2, vcut_left_c2], [0, size], c = 'deeppink')
            axs[1, 0].plot([0, size], [hcut_top_c2, hcut_top_c2], c = 'deeppink')
            axs[1, 0].plot([0, size], [hcut_bottom_c2, hcut_bottom_c2], c = 'deeppink')
            axs[0, 1].imshow(cam1[hcut_bottom_c1:hcut_top_c1, vcut_left_c1:vcut_right_c1], origin = 'lower', cmap = 'gray')
            axs[1, 1].imshow(cam2[hcut_bottom_c2:hcut_top_c2, vcut_left_c2:vcut_right_c2], origin = 'lower', cmap = 'gray')
            axs[0, 0].set_xlim(0, size)
            axs[0, 0].set_ylim(0, size)
            axs[1, 0].set_xlim(0, size)
            axs[1, 0].set_ylim(0, size)    
            axs[0, 0].set_ylabel("Cam 1")
            axs[1, 0].set_ylabel("Cam 2")

            plt.tight_layout()
            plt.show()

        if verbose:
            print(f"Cam 1 Field stop: {cam1_fieldstop}")
            print(f"Cam 2 Field stop: {cam2_fieldstop}")

        print(f"Field stop computation finished in {round(time.time() - tic, 3)}s.")
        return cam1_fieldstop, cam2_fieldstop

    else:
        raise Exception("Please choose 'fixed' or 'auto' fieldstop mmode")

# TESTING
"""
image_path_c1 = "/home/pablo/Desktop/SuObsTEsts/DATA/27/2024_05_01_10_10_36_149_0_690.img"
image_path_c2 = "/home/pablo/Desktop/SuObsTEsts/DATA/27/2024_05_01_10_10_36_149_1_690.img"
c1, _ = read_Tumag(image_path_c1)
c2, _ = read_Tumag(image_path_c2)

c1_fs, c2_fs = find_fieldstop(mode = "auto", cam1 = c1, cam2 = c2, plot_flag= True)
"""

def fieldstopping_and_shifting(cam1, cam2, fsc1, fsc2):

    new_cam1 = np.zeros((np.shape(cam1)))
    new_cam2 = np.zeros((np.shape(cam1)))

    field_stop_mask_cam1 = np.zeros((np.shape(cam1)))
    field_stop_mask_cam2 = np.zeros((np.shape(cam2)))

    field_stop_mask_cam1[fsc1[0][0] : fsc1[0][1], fsc1[1][0] : fsc1[1][1]] = 1

    new_cam1 = field_stop_mask_cam1 * cam1

    new_cam2[fsc1[0][0] : fsc1[0][1], fsc1[1][0] : fsc1[1][1]] = cam2[fsc2[0][0] : fsc2[0][1], fsc2[1][0] : fsc2[1][1]]    

    return new_cam1, new_cam2





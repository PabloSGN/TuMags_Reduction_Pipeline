"""
.. py:module:: TuMagTools.MicroPols_Extractor
.. module:: utils
        :platform: Unix
        :synopsis: Module to analyze micropols observations. Refer to jupyter 
        notebook : Long_cal_analysis.ipynb for example of use. 

        Command line execution : No command line execution.

.. moduleauthor:: Pablo Santamarina  (SPG - IAA) <psanta@iaa.es>
"""
# ============================ IMPORTS ====================================== #

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Own libs
from utils import read_Tumag
 
# ============================= CONFIG ====================================== #

plt.style.use('default')

# Index for each micropol in each of the 3 x 3 sets
MpolsConfig = {
                90 : [0, 5, 7],
                0 : [2, 3],
                45 : [4, 6], 
                135 : [1, 8] 
                }

Vlines = [] # Vertical division of the sets
Hlines = [] # Horizontal division of the sets
for i in range(16):
    Vlines.append(225 + 105 * i)
    Hlines.append(230 + 105 * i)

tilt_of_lines = 8 # tilt of the edges of the micropols 


Boxsize = 35 # Size of each MicroPol

# =========================================================================== #

def plot_square(x1, x2, y1, y2, axs, col):

    axs.plot([x1, x1], [y1, y2], ls = '--', c = col)
    axs.plot([x2, x2], [y1, y2], ls = '--', c = col)
    axs.plot([x1, x2], [y1, y1], ls = '--', c = col)
    axs.plot([x1, x2], [y2, y2], ls = '--', c = col)

def extract_modulations(box, boxsize):
    count = 0
    mpols = []
    for i in range(3):
        for j in range(3):    
            mpl = box[boxsize * j : boxsize * (j + 1), boxsize * i: boxsize * (i + 1)]
            mpols.append(mpl[5:30, 5:30])
            count += 1 

    return np.array(mpols)

def extract_micropols(image, Vlines, Hlines, tilt_px = tilt_of_lines):

    def tilt(x, tilt_px):
        return x * tilt_px/504 
    
    count = 0
    MP_Sets = {}
    for i in range(4):
        for j in range(4):

            x1 = Vlines[i] + tilt(Hlines[j], tilt_px)
            x2 = Vlines[i + 1] + tilt(Hlines[j + 1], tilt_px)
            y1 = Hlines[j] - tilt(Vlines[i], tilt_px)
            y2 = Hlines[j + 1] - tilt(Vlines[i + 1], tilt_px)

            MP_Sets[count] = {"set" : image[round(y1) : round(y2), round(x1) : round(x2)],
                              "mpols" : extract_modulations(image[round(y1) : round(y2), round(x1) : round(x2)], 35)}
            
            count += 1

    return MP_Sets

def get_mean_value(polangle, sets, MP_sets):

    Mean = []
    pols = []
    for st in sets:
        for index in MpolsConfig[polangle]:
            pols.append(MP_sets[st]["mpols"][index])

    pols = np.array(pols)

    return np.mean(pols)


def plot_configuration(image, selected_sets, polangle):

    fig, axs = plt.subplots(figsize = (12, 12))
    axs.imshow(image, cmap = 'inferno')

    for i in range(len(Vlines)):
        axs.plot([Vlines[i], Vlines[i] + (8 * 1800)/ 504], [Hlines[0], 1800], c = 'w', ls = '-', lw = 2)
        axs.plot([Vlines[0], 1800], [Hlines[i], Hlines[i] - (8 * 1800)/ 504], c = 'w', ls = '-', lw = 2)
        
        """axs.plot([Vlines[i] + Boxsize * 1, Vlines[i] + Boxsize * 1 + (8 * 1800)/ 504 ], [Hlines[0], 1800], c = 'w', ls = '--', lw = 1)
        axs.plot([Vlines[i] + Boxsize * 2, Vlines[i] + Boxsize * 2 + (8 * 1800)/ 504 ], [Hlines[0], 1800], c = 'w', ls = '--', lw = 1)

        axs.plot([Vlines[0], 1800], [Hlines[i] + Boxsize, Hlines[i] + Boxsize - (8 * 1800)/ 504], [Hlines[0], 1800], c = 'w', ls = '--', lw = 1)
        axs.plot([Vlines[0], 1800], [Hlines[i] + Boxsize * 2, Hlines[i] + Boxsize * 2 - (8 * 1800)/ 504], [Hlines[0], 1800], c = 'w', ls = '--', lw = 1)"""
    
    def tilt(x, tilt_px):
        return x * tilt_px/504 

    """coords = []
    tilt_px = tilt_of_lines
    count = 0
    for i in range(4):
        for j in range(4):
            x1 = Vlines[i] + tilt(Hlines[j], tilt_px)
            x2 = Vlines[i + 1] + tilt(Hlines[j + 1], tilt_px)
            y1 = Hlines[j] - tilt(Vlines[i], tilt_px)
            y2 = Hlines[j + 1] - tilt(Vlines[i + 1], tilt_px)

            axs.text(x1 + 45, y1 + 45, f"S{count}", color = 'w', fontsize = 'large')
            coords.append([x1, y1])
            count += 1

    for st in selected_sets:
        rect = patches.Rectangle((coords[st][0], coords[st][1]), 105, 105,
        facecolor='crimson', alpha = 0.5)

        count = 0
        boxsize = 35
        for i in range(3):
            for j in range(3):    
                if count in MpolsConfig[polangle]:
                    plot_square(boxsize * i + coords[st][0], boxsize * (i + 1) + coords[st][0],
                                 boxsize * j + coords[st][1], boxsize * (j + 1)+ coords[st][1], 
                                 axs, 'k')
                count += 1 


    # Add the patch to the Axes
        axs.add_patch(rect)"""

  

def get_modulation_curve(Folder, selected_sets, polangle):
    
    all_images = sorted(glob.glob(os.path.join(Folder, "*img")))
    print(f"Images found : {len(all_images)}")
    print(f"Selected Sets: {selected_sets}")
    print(f"Selected pol angle: {polangle}")

    mean = np.zeros((2, 25))
    volts = np.zeros(25)
    for img_ind, img in enumerate(all_images):
        I, H = read_Tumag(img)

        mpols = extract_micropols(I, Vlines, Hlines)

        if img_ind < 25:
            volts[img_ind] = H["Rocli1_LCVR"] * 0.00179439
            mean[0, img_ind] = get_mean_value(polangle, selected_sets, mpols)
        else:
            mean[1, img_ind - 25] = get_mean_value(polangle, selected_sets, mpols)

    fig, axs = plt.subplots(figsize = (15, 8))
    axs.plot(volts, mean[0], c = 'deeppink', lw = 3, label = 'I + V')
    axs.plot(volts, mean[1], c = 'gold', lw = 3, label = 'I - V')
    axs.set_xlabel("LCVR1 Volts")
    axs.set_ylabel("Intensity")
    axs.set_title(f"Polarization angle: {polangle}")
    axs.grid(True, c = 'w', alpha = 0.2)
    axs.legend()








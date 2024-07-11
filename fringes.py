# ---------------------------- DESCRIPTION --------------------------------------- #

# ------------------------------ IMPORTS ----------------------------------------- #

import os
import glob 
import numpy as np
import matplotlib.pyplot as plt

from field_stop_finder import find_fieldstop

# ------------------------------ CONFIG ------------------------------------------ #

# ------------------------------  CODE  ------------------------------------------ # 

def fringe_freq_finder(img):

    fft = np.fft.fft2(img)

    freqs = np.fft.fftshift(np.abs(np.fft.fft2(img)))
    fig, axs = plt.subplots(1, 2, figsize = (10, 5))

    axs[0].imshow(img, cmap = 'inferno')
    axs[1].imshow(freqs, cmap = 'afmhot', norm = 'log')

    plt.show()

# img_example = np.load("test_ff.npy")
# fs_c1, fs_c2 = find_fieldstop("fixed")
# fringe_freq_finder(img_example[1, fs_c2[0][0]:fs_c2[0][1], fs_c2[1][0]:fs_c2[1][1]])






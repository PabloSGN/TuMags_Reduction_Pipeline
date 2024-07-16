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
import os 
import glob

from scipy.signal import correlate, correlate2d

import matplotlib.pyplot as plt

# Own libs
from utils import read_Tumag
import config as cf
from field_stop_finder import find_fieldstop, fieldstopping_and_shifting

# ------------------------------  CODE  ------------------------------------------ # 


pinholes1 = "/home/pablo/Desktop/TuMAGDATA/flare_set_ordered/Spectral_calibration/36"
pinholes2 = "/home/pablo/Desktop/TuMAGDATA/flare_set_ordered/Spectral_calibration/37"
pinholes3 = "/home/pablo/Desktop/TuMAGDATA/flare_set_ordered/Spectral_calibration/38"
ims_c1 = sorted(glob.glob(f"{pinholes1}/*_0_*"))
ims_c2 = sorted(glob.glob(f"{pinholes1}/*_1_*"))

flats = "/home/pablo/Desktop/TuMAGDATA/flare_set_ordered/2.02/39"
flats_c1 = sorted(glob.glob(f"{flats}/*_0_*"))
flats_c2 = sorted(glob.glob(f"{flats}/*_1_*"))

i1, _ = read_Tumag(ims_c1[0])
i2, _ = read_Tumag(ims_c2[0])

i1 = i1 / np.max(i1)
i2 = i2 / np.max(i2)
i2_flipped = np.flip(i2, axis = -1)

f1, _ = read_Tumag(flats_c1[0])
f2, _ = read_Tumag(flats_c2[0])

f1 = f1 / np.max(f1)
f2 = f2 / np.max(f2)
f2_flipped = np.flip(f2, axis = -1)

#corr_ph = correlate(i1, i2_flipped, mode='full', method = 'fft')
#corr_flats = correlate(f1, f2_flipped, mode='full', method = 'fft')

stopx = [875, 1025]
stopy = [925, 1075]

corr_ph = correlate(i1[stopx[0]:stopx[1], stopy[0]:stopy[1]], 
                    i2_flipped[stopx[0]:stopx[1], stopy[0]:stopy[1]], mode='full')

max_index = np.unravel_index(np.argmax(corr_ph), corr_ph.shape)
shift_x = max_index[0] - (stopx[1] - stopx[0] - 1)
shift_y = max_index[1] - (stopy[1] - stopy[0] - 1)
print("shift : ", shift_x, shift_y)

fig, axs = plt.subplots(1, 4, figsize = (15, 5))

axs[0].imshow(i1[stopx[0]:stopx[1], stopy[0]:stopy[1]], cmap = "afmhot")
axs[1].imshow(i2_flipped[stopx[0]:stopx[1], stopy[0]:stopy[1]], cmap = "afmhot")
im = axs[2].imshow(i1[stopx[0]:stopx[1], stopy[0]:stopy[1]]- i2_flipped[stopx[0]:stopx[1], stopy[0]:stopy[1]], cmap = "inferno")
plt.colorbar(im,fraction=0.046, pad=0.04)

fs_c1  = find_fieldstop(cam1 = f1, plot_flag = False, verbose = True)


fs_c2 = np.zeros(np.shape(fs_c1)).astype("int")

fs_c2[0] = fs_c1[0] - int(shift_x)
fs_c2[1] = fs_c1[1] - int(shift_y)

print(fs_c1, fs_c2)

ph_shift1, ph_shift2 = fieldstopping_and_shifting(i1, i2_flipped, fs_c1, fs_c2)

im = axs[3].imshow(ph_shift1[stopx[0]:stopx[1], stopy[0]:stopy[1]]- ph_shift2[stopx[0]:stopx[1], stopy[0]:stopy[1]], cmap = "inferno")
plt.colorbar(im,fraction=0.046, pad=0.04)

plt.tight_layout()
plt.show()

"""
# Find the peak in the cross-correlation
max_index = np.unravel_index(np.argmax(corr_flats), corr_flats.shape)

# Calculate the shifts
shift_y = max_index[0] - (2016 - 1)
shift_x = max_index[1] - (2016 - 1)
print(shift_x, shift_y)

# corr_flat = correlate(f1, f2_flipped, mode = 'same', method = 'fft')
# corr_phol = correlate(i1, i2_flipped, mode = 'same', method = 'fft')


fig, axs = plt.subplots(2, 3, figsize = (10, 5))

axs[0, 0].set_ylabel("Pre_alignment")
axs[1, 0].set_ylabel("Pre_alignment")


axs[0, 0].set_title("Original - shift")
im = axs[0, 0].imshow(i2_flipped - i2, cmap = "afmhot", vmin = 0, vmax = 0.05)
plt.colorbar(im,fraction=0.046, pad=0.04)

axs[0, 1].set_title("Ph Diff")
im = axs[0, 1].imshow(i1 - i2_flipped, cmap = "afmhot", vmin = 0, vmax = 0.05)
plt.colorbar(im,fraction=0.046, pad=0.04)

axs[0, 2].set_title("Flat diff")
im = axs[0, 2].imshow(f1*0.5 - f2_flipped*0.5, cmap = "afmhot")#, vmin = 2000, vmax = 2400)
plt.colorbar(im,fraction=0.046, pad=0.04)

fs_c1, _ = find_fieldstop(mode = "auto", cam1 = f1, cam2 = f2_flipped, plot_flag = False, verbose = True)


fs_c2= [[fs_c1[0][0] - shift_y, fs_c1[0][1] -shift_y], 
        [fs_c1[1][0] - shift_x, fs_c1[1][1] - shift_x]]

f1_shift, f2_shift = fieldstopping_and_shifting(f1, f2_flipped, fs_c1, fs_c2)
ph1_shift, ph2_shift = fieldstopping_and_shifting(i1, i2_flipped, fs_c1, fs_c2)


im = axs[1, 0].imshow(f2 - f2_shift, cmap = "afmhot", vmin = 0, vmax = 0.2)
plt.colorbar(im,fraction=0.046, pad=0.04)

im = axs[1, 1].imshow(ph1_shift - ph2_shift, cmap = "afmhot", vmin = 0, vmax = 0.05)
plt.colorbar(im,fraction=0.046, pad=0.04)

im = axs[1, 2].imshow(f1_shift - f2_shift, cmap = "afmhot")#, vmin = 2000, vmax = 2400)
plt.colorbar(im,fraction=0.046, pad=0.04)

plt.tight_layout()
plt.show()"""





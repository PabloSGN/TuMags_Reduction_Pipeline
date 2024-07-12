# ---------------------------- DESCRIPTION --------------------------------------- #

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in 
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# Own Libs
from utils import read_Tumag
from field_stop_finder import find_fieldstop
from master_dark import compute_master_darks
from master_flatfield import compute_master_flat_field

# ------------------------------ CONFIG ------------------------------------------ #

# Testeo
daks_c1 = sorted(glob.glob(f"/home/pablo/Desktop/TuMAGDATA/Darks/*_0_*")) 
daks_c2 = sorted(glob.glob(f"/home/pablo/Desktop/TuMAGDATA/Darks/*_1_*")) 

image_path_c1 = "/home/pablo/Desktop/SuObsTEsts/DATA/27/2024_05_01_10_10_36_149_0_690.img"
image_path_c2 = "/home/pablo/Desktop/SuObsTEsts/DATA/27/2024_05_01_10_10_36_149_1_690.img"

flats_c1 = sorted(glob.glob(f"/home/pablo/Desktop/SuObsTEsts/DATA/27/*_0_*")) 
flats_c2 = sorted(glob.glob(f"/home/pablo/Desktop/SuObsTEsts/DATA/27/*_1_*")) 

# ------------------------------  CODE  ------------------------------------------ # 

c1, _ = read_Tumag(image_path_c1)
c2, _ = read_Tumag(image_path_c2)

fs_c1, fs_c2 = find_fieldstop(mode = "auto", cam1 = c1, cam2 = c2, plot_flag = False, verbose = True)

dc = compute_master_darks(daks_c1, daks_c2, fs_c1, fs_c2, verbose = True)

ff = compute_master_flat_field(flats_c1, flats_c2, dc, fs_c1, fs_c2, verbose = True)

fig, axs = plt.subplots(2, 2, figsize = (10, 10))

axs[0, 0].imshow(dc[0], cmap = "gray")
axs[0, 1].imshow(dc[1], cmap = "gray")

axs[1, 0].imshow(ff[0], cmap = "inferno")
axs[1, 1].imshow(ff[1], cmap = "inferno")

np.save("test_ff.npy", ff)

plt.tight_layout()
plt.show()
    

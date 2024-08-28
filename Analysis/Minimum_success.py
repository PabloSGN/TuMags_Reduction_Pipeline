# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in 
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Own Libs
sys.path.append("/home/users/dss/orozco/Tumag/PabloTests")
import config as cf
from utils import read_Tumag
from field_stop_finder import compute_alignment, apply_fieldstop_and_align_array
from master_dark import compute_master_darks
from master_flatfield import compute_master_flat_field
import image_handler as ih
from demodulation import demodulate

# ------------------------------ CONFIG -------------------------------------------- # 

# Indexes. 
flats_obs1_index = "D10-2740-4339"
flats_obs2_index = "D10-2740-4339"
darks_index = "D10-5620-5719"
pinholes_index = "D10-6188-6193"
obs_index = "D10-304-2439"

dark_paths = ih.get_images_paths(darks_index)
ff_obs1_paths = ih.get_images_paths(flats_obs1_index)
pinholes_paths = ih.get_images_paths(pinholes_index)
obs_images = ih.get_images_paths(obs_index)

# ----------------------------------------------------------------------------------- #

# Compute dark current

dc = compute_master_darks(dark_paths[:-6], verbose = True) # There is something strange in the last darks, you can see some structure -> Me los cargo con el :-6


# Compute field-stop and alignment.
# We need flat_fields for field-stop and Pinholes for alignment. 
flatc1, h = read_Tumag([x for x in ff_obs1_paths if "_0_" in x][0])
flatc2, _ = read_Tumag([x for x in ff_obs1_paths if "_1_" in x][0])

print(f"\n Nacc of flats: {h['nAcc']}")

fs_c1, fs_c2 = compute_alignment(flat_cam1=flatc1 / np.max(flatc1), flat_cam2= flatc2 / np.max(flatc2), pinhole_c1_path=pinholes_paths[0], pinhole_c2_path=pinholes_paths[1], method = 'pinhole', plot_flag=False, verbose = True)

# Compute master flat field
ff = compute_master_flat_field(ff_obs1_paths, dc = dc, verbose = True)[0]

# Align flats
ff = apply_fieldstop_and_align_array(ff, fs_c1, fs_c2)

# Images of first mode
first_mode = obs_images[:cf.om_config["1"]["images_per_mode"]]

# Read first image to get accumulations
_, H = ih.read(first_mode[0])

print("\nProcessing images of observing mode...")
ob = ih.nominal_observation("1", first_mode, dc * H["nAcc"])

# Align observation
data = apply_fieldstop_and_align_array(ob.get_data(), fs_c1, fs_c2)
info = ob.get_info()

# Correct flats
corrected = np.zeros(np.shape(data))
for mod in range(info["Nmods"]):
    for lamb in range(info["Nlambda"]):

        corrected[1, lamb, mod] = (data[0, lamb, mod]) / ff[0, lamb, mod]
        corrected[0, lamb, mod] = (data[1, lamb, mod]) / ff[1, lamb, mod]


fig, axs = plt.subplots(2, 4, figsize = (20, 10))

im = axs[0, 0].imshow(dc[0], cmap = 'gray', vmin = np.mean(dc) - np.std(dc), vmax = np.mean(dc) + np.std(dc))
divider = make_axes_locatable(axs[0, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[1, 0].imshow(dc[1], cmap = 'gray', vmin = np.mean(dc) - np.std(dc), vmax = np.mean(dc) + np.std(dc))
divider = make_axes_locatable(axs[1, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)


im = axs[0, 1].imshow(ff[0,-1, 0], cmap = 'inferno')
divider = make_axes_locatable(axs[0, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[1, 1].imshow(ff[0,-1, 0], cmap = 'inferno')
divider = make_axes_locatable(axs[1, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

im = axs[0, 2].imshow(corrected[0,-1, 0,  fs_c1[0][0]:fs_c1[0][1], fs_c1[1][0]:fs_c1[1][1]], cmap = 'inferno')
divider = make_axes_locatable(axs[0, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[1, 2].imshow(corrected[0,-1, 0, fs_c1[0][0]:fs_c1[0][1], fs_c1[1][0]:fs_c1[1][1]], cmap = 'inferno')
divider = make_axes_locatable(axs[1, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

cds_check_x = [950 + 125, 1150] 
cds_check_y = [200, 400 - 125] 

im = axs[0, 3].imshow(corrected[0, -1,0, cds_check_x[0]:cds_check_x[1], cds_check_y[0]:cds_check_y[1]] / np.max(corrected[0, -1,0, cds_check_x[0]:cds_check_x[1], cds_check_y[0]:cds_check_y[1]]) -\
                   corrected[1, -1,0, cds_check_x[0]:cds_check_x[1], cds_check_y[0]:cds_check_y[1]] / np.max(corrected[1, -1,0, cds_check_x[0]:cds_check_x[1], cds_check_y[0]:cds_check_y[1]]), cmap = "Spectral")

divider = make_axes_locatable(axs[0, 3])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.show()

# ------------------------------ IMPORTS ----------------------------------------- #

import sys
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
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
import phase_diversity as pd

# ------------------------------ CONFIG -------------------------------------------- # 

pd_indexes = "D10-6194-6513"
flats_mode1_index = "D10-2740-4339"
flats_mode2_index = "D10-4340-5619"
dark_indexes = "D10-5620-5719"
pinhole_indexes = "D10-6188-6193"

# Get the path to the images
dark_paths = ih.get_images_paths(dark_indexes)
ff_obs1_paths = ih.get_images_paths(flats_mode1_index)
ff_obs2_paths = ih.get_images_paths(flats_mode2_index)
pinholes_paths = ih.get_images_paths(pinhole_indexes)
# PD has its own function to parse images. 

# ----------------------------------------------------------------------------------- #

# Start by computing master dark
dc = compute_master_darks(dark_paths[:-6], verbose = True) # There is something strange in the last darks, you can see some structure -> Me los cargo con el :-6

# Compute flat-fields
# Compute master flat field
ff_obs1, ff_obs1_info = compute_master_flat_field(ff_obs1_paths, dc = dc, verbose = True)
ff_obs2, ff_obs2_info = compute_master_flat_field(ff_obs2_paths, dc = dc, verbose = True)

# Select continuum wavelength and modulation 0 and normalize
ff_1 = ff_obs1[:, -1, 0] / np.max(ff_obs1[:, -1, 0])
ff_2 = ff_obs2[:, -1, 0] / np.max(ff_obs2[:, -1, 0])

# flat field both cams, cont wave, modul 0. 
pd_data_fe = pd.process_pd_observation(pd_indexes, filt = 0, verbose = True) 
pd_data_mg = pd.process_pd_observation(pd_indexes, filt = 1, verbose = True) 


print("Plotting...")

fig, axs = plt.subplots(2, 3, figsize = (10, 10))
im = axs[0, 0].imshow(ff_1[0], cmap = "magma")
divider = make_axes_locatable(axs[0, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[0, 1].imshow(ff_2[0], cmap = "magma")
divider = make_axes_locatable(axs[0, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[1, 0].imshow(ff_1[1], cmap = "magma")
divider = make_axes_locatable(axs[1, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[1, 1].imshow(ff_2[1], cmap = "magma")
divider = make_axes_locatable(axs[1, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[0, 2].imshow(dc[0], cmap = "magma")
divider = make_axes_locatable(axs[0, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[1, 2].imshow(dc[1], cmap = "magma")
divider = make_axes_locatable(axs[1, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

plt.tight_layout()
plt.savefig("flats_n_darks.png", bbox_inches = 'tight')



fig, axs = plt.subplots(2, 4, figsize = (17, 10))
im = axs[0, 0].set_title("Iron defocused")
im = axs[0, 1].set_title("Iron focused")
im = axs[0, 2].set_title("Magnesium defocused")
im = axs[0, 3].set_title("Magnesium focused")

im = axs[0, 0].imshow(pd_data_fe[0, 0, 1], cmap = "gray", vmin = np.mean(pd_data_fe[0, 0, 1]) - np.std(pd_data_fe[0, 0, 1]) * 2,
                                                          vmax = np.mean(pd_data_fe[0, 0, 1]) + np.std(pd_data_fe[0, 0, 1]) * 2 )
divider = make_axes_locatable(axs[0, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[0, 1].imshow(pd_data_fe[0, 1, 1], cmap = "gray", vmin = np.mean(pd_data_fe[0, 1, 0]) - np.std(pd_data_fe[0, 1, 0]) * 2,
                                                          vmax = np.mean(pd_data_fe[0, 1, 0]) + np.std(pd_data_fe[0, 1, 0]) * 2 )
divider = make_axes_locatable(axs[0, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[0, 2].imshow(pd_data_mg[0, 0, 1], cmap = "gray", vmin = np.mean(pd_data_mg[0, 0, 0]) - np.std(pd_data_mg[0, 0, 0]) * 2,
                                                          vmax = np.mean(pd_data_mg[0, 0, 0]) + np.std(pd_data_mg[0, 0, 0]) * 2 )
divider = make_axes_locatable(axs[0, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[0, 3].imshow(pd_data_mg[0, 1, 1], cmap = "gray", vmin = np.mean(pd_data_mg[0, 1, 0]) - np.std(pd_data_mg[0, 1, 0]) * 2,
                                                          vmax = np.mean(pd_data_mg[0, 1, 0]) + np.std(pd_data_mg[0, 1, 0]) * 2 )
divider = make_axes_locatable(axs[0, 3])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[1, 0].imshow(pd_data_fe[1, 0, 1], cmap = "gray", vmin = np.mean(pd_data_fe[1, 0, 0]) - np.std(pd_data_fe[1, 0, 0]) * 2,
                                                          vmax = np.mean(pd_data_fe[1, 0, 0]) + np.std(pd_data_fe[1, 0, 0]) * 2 )
divider = make_axes_locatable(axs[1, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[1, 1].imshow(pd_data_fe[1, 1, 1], cmap = "gray", vmin = np.mean(pd_data_fe[1, 1, 0]) - np.std(pd_data_fe[1, 1, 0]) * 2,
                                                          vmax = np.mean(pd_data_fe[1, 1, 0]) + np.std(pd_data_fe[1, 1, 0]) * 2 )
divider = make_axes_locatable(axs[1, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[1, 2].imshow(pd_data_mg[1, 0, 1], cmap = "gray", vmin = np.mean(pd_data_mg[1, 0, 0]) - np.std(pd_data_mg[1, 0, 0]) * 2,
                                                          vmax = np.mean(pd_data_mg[1, 0, 0]) + np.std(pd_data_mg[1, 0, 0]) * 2 )
divider = make_axes_locatable(axs[1, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[1, 3].imshow(pd_data_mg[1, 1, 1], cmap = "gray", vmin = np.mean(pd_data_mg[1, 1, 0]) - np.std(pd_data_mg[1, 1, 0]) * 2,
                                                          vmax = np.mean(pd_data_mg[1, 1, 0]) + np.std(pd_data_mg[1, 1, 0]) * 2 )
divider = make_axes_locatable(axs[1, 3])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)


im = axs[0, 0].set_ylabel("CAM 1")
im = axs[1, 0].set_ylabel("CAM 2")
plt.tight_layout()
plt.savefig("PD_example.png", bbox_inches = 'tight')
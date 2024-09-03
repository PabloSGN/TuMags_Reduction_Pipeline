# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in 
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle

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
flats_obs2_index = "D10-4340-5619"
darks_index = "D10-5620-5719"
pinholes_index = "D10-6188-6193"
obs_index = "D10-304-2439"

dark_paths = ih.get_images_paths(darks_index)
ff_obs1_paths = ih.get_images_paths(flats_obs1_index)
ff_obs2_paths = ih.get_images_paths(flats_obs2_index)
pinholes_paths = ih.get_images_paths(pinholes_index)
obs_images = ih.get_images_paths(obs_index)

# ----------------------------------------------------------------------------------- #

# Compute dark current
dc = compute_master_darks(dark_paths[:-6], verbose = True) # There is something strange in the last darks, you can see some structure -> Me los cargo con el :-6

# Plotting Dark current
fig, axs = plt.subplots(1, 2, figsize = (14, 7))
im1 = axs[0].imshow(dc[0], cmap = 'gray', vmin = np.mean(dc[0]) - np.std(dc[0]) * 2, vmax = np.mean(dc[0]) + np.std(dc[0]) * 2)
im2 = axs[1].imshow(dc[1], cmap = 'gray', vmin = np.mean(dc[0]) - np.std(dc[0]) * 2, vmax = np.mean(dc[0]) + np.std(dc[0]) * 2)

axs[0].set_title("Cam 1")
axs[0].set_title("Cam 2")

# Colorbars
divider = make_axes_locatable(axs[0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im1, cax=cax)

divider = make_axes_locatable(axs[1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im2, cax=cax)
plt.tight_layout()

plt.savefig(f"/home/users/dss/orozco/Tumag/PabloTests/Analysis/Minimum_success_plots/DarkCurrent.png", bbox_inches = "tight")

# ----------------------------------------------------------------------------------- #
# Compute flat-fields


ff_obs1, ff_obs1_info = compute_master_flat_field(ff_obs1_paths, dc = dc, verbose = True)
ff_obs2, ff_obs2_info = compute_master_flat_field(ff_obs2_paths, dc = dc, verbose = True)

# Compute field-stop and alignment.

fs_c1, fs_c2 = compute_alignment(flat_cam1 = ff_obs1[0, -1, 0] / np.max(ff_obs1[0, -1, 0]),
                                 flat_cam2 = ff_obs1[1, -1, 0] / np.max(ff_obs1[1, -1, 0]),
                                 pinhole_c1_path=pinholes_paths[0], pinhole_c2_path=pinholes_paths[1], method = 'pinhole', plot_flag=False, verbose = True)

# Align flats
print("Aligning flats")
ff1 = apply_fieldstop_and_align_array(ff_obs1, fs_c1, fs_c2)
ff2 = apply_fieldstop_and_align_array(ff_obs2, fs_c1, fs_c2)

# ----------------------------------------------------------------------------------- #

"""
Observation_counters = ih.separate_ocs(obs_images, verbose = True)
with open("minimum_success_observation_counters.pickle", 'wb') as file:
    # Serialize and write the variable to the file
    pickle.dump(Observation_counters, file)
"""
    
# Open the file
with open("minimum_success_observation_counters.pickle", 'rb') as file:
    Observation_counters = pickle.load(file)

print("\nProcessing images of observing mode...")
ob_mg = ih.nominal_observation("1", Observation_counters[81], dc)
print("Mode 1 processed.")

ob_fe = ih.nominal_observation("2.02", Observation_counters[82], dc)
print("Mode 2.02 processed.")


print("Plotting")
fig, axs = plt.subplots(1, 2, figsize = (14, 7))
axs[0].plot(cf.om_config["1"]["lambda_array"], np.mean(ob_mg.get_data()[0, :, 0], axis = (1, 2)), marker = 'x', c = 'indigo', lw = 2, label = "Mod. 1")
axs[0].plot(cf.om_config["1"]["lambda_array"], np.mean(ob_mg.get_data()[0, :, 1], axis = (1, 2)), marker = 'x', c = 'darkorange', lw = 2, label = "Mod. 2")
axs[0].plot(cf.om_config["1"]["lambda_array"], np.mean(ob_mg.get_data()[0, :, 2], axis = (1, 2)), marker = 'x', c = 'crimson', lw = 2, label = "Mod. 3")
axs[0].plot(cf.om_config["1"]["lambda_array"], np.mean(ob_mg.get_data()[0, :, 3], axis = (1, 2)), marker = 'x', c = 'mediumseagreen', lw = 2, label = "Mod. 4")
axs[1].plot(cf.om_config["2.02"]["lambda_array"], np.mean(ob_fe.get_data()[0, :, 0], axis = (1, 2)), marker = 'x', c = 'indigo', lw = 2)
axs[1].plot(cf.om_config["2.02"]["lambda_array"], np.mean(ob_fe.get_data()[0, :, 1], axis = (1, 2)), marker = 'x', c = 'darkorange', lw = 2)
axs[1].plot(cf.om_config["2.02"]["lambda_array"], np.mean(ob_fe.get_data()[0, :, 2], axis = (1, 2)), marker = 'x', c = 'crimson', lw = 2)
axs[1].plot(cf.om_config["2.02"]["lambda_array"], np.mean(ob_fe.get_data()[0, :, 3], axis = (1, 2)), marker = 'x', c = 'mediumseagreen', lw = 2)
axs[0].legend(edgecolor = 'k')
axs[1].set_xlabel(r"$\Delta \lambda$ [m$\AA$]")
axs[1].set_ylabel("Average intensity")
axs[0].set_xlabel(r"$\Delta \lambda$ [m$\AA$]")
axs[0].set_ylabel("Average intensity")
axs[0].set_title("Mode 1")
axs[1].set_title("Mode 2.02")
axs[0].grid(True, c = 'k', alpha = 0.3)
axs[1].grid(True, c = 'k', alpha = 0.3)
plt.tight_layout()
plt.savefig(f"/home/users/dss/orozco/Tumag/PabloTests/Analysis/Minimum_success_plots/Average_profile.png", bbox_inches = "tight")

print("Finished first plot")

# Align observation
data_mg = apply_fieldstop_and_align_array(ob_mg.get_data(), fs_c1, fs_c2)
info_mg = ob_mg.get_info()

mask = ff1[0, 0, 0] != 0
# Correct flats
corrected_mg = np.zeros(np.shape(data_mg))
for mod in range(info_mg["Nmods"]):
    for lamb in range(info_mg["Nlambda"]):
        corrected_mg[1, lamb, mod, mask] = (data_mg[0, lamb, mod, mask]) / ff1[0, lamb, mod, mask]
        corrected_mg[0, lamb, mod, mask] = (data_mg[1, lamb, mod, mask]) / ff1[1, lamb, mod, mask]


# Align observation
data_fe = apply_fieldstop_and_align_array(ob_fe.get_data(), fs_c1, fs_c2)
info_fe = ob_fe.get_info()

# Correct flats
corrected_fe = np.zeros(np.shape(data_fe))
for mod in range(info_fe["Nmods"]):
    for lamb in range(info_fe["Nlambda"]):
        corrected_fe[1, lamb, mod, mask] = (data_fe[0, lamb, mod, mask]) / ff2[0, lamb, mod, mask]
        corrected_fe[0, lamb, mod, mask] = (data_fe[1, lamb, mod, mask]) / ff2[1, lamb, mod, mask]

        
fig, axs = plt.subplots(2, 2, figsize = (12, 12))
axs[0, 0].set_title("Magnesium flats")
im = axs[0, 0].imshow(ff1[0, -1, 0], cmap = 'inferno', vmin = np.mean(ff1[0, -1, 0]) - np.std(ff1[0, -1, 0]) * 3, vmax = np.mean(ff1[0, -1, 0]) + np.std(ff1[0, -1, 0]) * 3)
divider = make_axes_locatable(axs[0, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[1, 0].set_title("Iron flats")
im = axs[0, 1].imshow(ff2[0, -1, 0], cmap = 'inferno', vmin = np.mean(ff2[0, -1, 0]) - np.std(ff2[0, -1, 0]) * 3, vmax = np.mean(ff2[0, -1, 0]) + np.std(ff2[0, -1, 0]) * 3)
divider = make_axes_locatable(axs[0, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[0, 1].set_title("Mag. cont")
im = axs[0, 1].imshow(corrected_mg[0,-1, 0], cmap = 'gray')
divider = make_axes_locatable(axs[0, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[1, 1].set_title("Iron Cont")
axs[1, 1].imshow(corrected_fe[0,-1, 0], cmap = 'gray')
divider = make_axes_locatable(axs[1, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

plt.tight_layout()
plt.savefig(f"/home/users/dss/orozco/Tumag/PabloTests/Analysis/Minimum_success_plots/Flat_correction.png", bbox_inches = "tight")


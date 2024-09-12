# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in 
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from astropy.io import fits
import matplotlib as mp

# Own Libs
sys.path.append("../")
import config as cf
from utils import read_Tumag
from field_stop_finder import compute_alignment, apply_fieldstop_and_align_array
from master_dark import compute_master_darks
from master_flatfield import compute_master_flat_field
import image_handler as ih
from demodulation import demodulate

# ------------------------------ CONFIG -------------------------------------------- # 

# CONFIG
plt.style.use("default")
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"
 

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

# Dark current. 

dark_mean_c1 = []
dark_mean_c2 = []

for image in dark_paths:
    I, H = ih.read(image)
    
    if H["cam"] == 0:
        dark_mean_c1.append(np.mean(I[300:-300, 300:-300]))
    else:
        dark_mean_c2.append(np.mean(I[300:-300, 300:-300]))

fig, axs = plt.subplots(1, 2, figsize = (10.5, 5))

axs[0].plot(dark_mean_c1, color = 'magenta', lw = 3, marker = 'x', label = "Cam 1")
axs[1].plot(dark_mean_c1, color = 'darkorange', lw = 3, marker = 'x', label = "Cam 2")
axs[0].grid(True, c = 'k', alpha = 0.3)
axs[1].grid(True, c = 'k', alpha = 0.3)

axs[0].set_ylabel("Average counts")
axs[0].set_xlabel("Nº Image")
axs[1].set_ylabel("Average counts")
axs[1].set_xlabel("Nº Image")

axs[0].legend(edgecolor = 'k')
axs[1].legend(edgecolor = 'k')

plt.savefig(f"ms_presentation/dark_mean.png", bbox_inches = "tight")

fig, axs = plt.subplots(figsize = (5.5, 5))
im = axs.imshow(I, vmin = np.mean(I) - np.std(I) * 2, vmax = np.mean(I) + np.std(I) * 2, cmap = "gray")
divider = make_axes_locatable(axs)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs.set_xticks([])
axs.set_yticks([])

plt.tight_layout()
plt.savefig(f"ms_presentation/wrong_dark.png", bbox_inches = "tight")

dc = compute_master_darks(dark_paths[:-6], verbose = True) # There is something strange in the last darks, you can see some structure -> Me los cargo con el :-6

fig, axs = plt.subplots(1, 2, figsize = (14, 7))
im1 = axs[0].imshow(dc[0], cmap = 'gray', vmin = np.mean(dc[0]) - np.std(dc[0]) * 2, vmax = np.mean(dc[0]) + np.std(dc[0]) * 2)
im2 = axs[1].imshow(dc[1], cmap = 'gray', vmin = np.mean(dc[0]) - np.std(dc[0]) * 2, vmax = np.mean(dc[0]) + np.std(dc[0]) * 2)

axs[0].set_title("Cam 1")
axs[1].set_title("Cam 2")

axs[0].set_xticks([])
axs[0].set_yticks([])
axs[1].set_xticks([])
axs[1].set_yticks([])

# Colorbars
divider = make_axes_locatable(axs[0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im1, cax=cax)

divider = make_axes_locatable(axs[1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im2, cax=cax)
plt.tight_layout()
plt.savefig(f"ms_presentation/DarkCurrent.png", bbox_inches = "tight")

"""
ff_obs1, ff_obs1_info = compute_master_flat_field(ff_obs1_paths, dc = dc, verbose = True)
ff_obs2, ff_obs2_info = compute_master_flat_field(ff_obs2_paths, dc = dc, verbose = True)

# Compute field-stop and alignment.

fs_c1, fs_c2 = compute_alignment(flat_cam1 = ff_obs2[0, -1, 0] / np.max(ff_obs2[0, -1, 0]),
                                 flat_cam2 = ff_obs2[1, -1, 0] / np.max(ff_obs2[1, -1, 0]),
                                 pinhole_c1_path=pinholes_paths[0], pinhole_c2_path=pinholes_paths[1], method = 'pinhole', plot_flag=False, verbose = True)

np.save("field_stop.npy", np.array([fs_c1, fs_c2]))

# Align flats
print("Aligning flats")
ff1 = apply_fieldstop_and_align_array(ff_obs1, fs_c1, fs_c2)
ff2 = apply_fieldstop_and_align_array(ff_obs2, fs_c1, fs_c2)


hdu = fits.PrimaryHDU(ff1)
hdu.writeto("flats_mode_1_ms.fits")

hdu = fits.PrimaryHDU(ff2)
hdu.writeto("flats_mode_2_ms.fits")
"""

ff1 = fits.getdata("flats_mode_1_ms.fits")
ff2 = fits.getdata("flats_mode_2_ms.fits")

fs = np.load("field_stop.npy")
fs_c1 = fs[0]
fs_c2 = fs[1]


fig, axs = plt.subplots(1, 4, figsize = (20, 5))

axs[0].set_title("Flat cam 1 - OM : 1 \n Mod : 0 - Continuum")
im = axs[0].imshow(ff1[0, -1, 0], cmap = 'inferno')
divider = make_axes_locatable(axs[0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[1].set_title("Flat cam 2 - OM : 1 \n Mod : 0 - Continuum")
im = axs[1].imshow(ff1[1, -1, 0], cmap = 'inferno')
divider = make_axes_locatable(axs[1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[2].set_title("Flat cam 1 - OM : 2.02 \n Mod : 0 - Continuum")
im = axs[2].imshow(ff2[0, -1, 0], cmap = 'inferno')
divider = make_axes_locatable(axs[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[3].set_title("Flat cam 2 - OM : 2.02 \n Mod : 0 - Continuum")
im = axs[3].imshow(ff1[1, -1, 0], cmap = 'inferno')
divider = make_axes_locatable(axs[3])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)


for i in range(4):
    axs[i].set_xticks([])
    axs[i].set_yticks([])

plt.tight_layout()
plt.savefig(f"ms_presentation/flat-fields.png", bbox_inches = "tight")

"""
OCs = ih.separate_ocs(obs_images, verbose = True)

with open("ms_ocs.pickle", 'wb') as file:
    # Serialize and write the variable to the file
    pickle.dump(OCs, file)"""

# Open the file
with open("ms_ocs.pickle", 'rb') as file:
    OCs = pickle.load(file)

"""
ob_mg = ih.nominal_observation("1", OCs[83]["ims"], dc)
ob_fe = ih.nominal_observation("2.02", OCs[82]["ims"], dc)

# Align observation
data_fe = apply_fieldstop_and_align_array(ob_fe.get_data(), fs_c1, fs_c2)
info_fe = ob_fe.get_info()

# Align observation
data_mg = apply_fieldstop_and_align_array(ob_mg.get_data(), fs_c1, fs_c2)
info_mg = ob_mg.get_info()

mask = ff2[0, 0, 0] != 0

# Correct flats
corrected_mg = np.zeros(np.shape(data_mg))
for mod in range(info_mg["Nmods"]):
    for lamb in range(info_mg["Nlambda"]):
        corrected_mg[0, lamb, mod, mask] = (data_mg[0, lamb, mod, mask]) / ff1[0, lamb, mod, mask]
        corrected_mg[1, lamb, mod, mask] = (data_mg[1, lamb, mod, mask]) / ff1[1, lamb, mod, mask]

# Correct flats
corrected_fe = np.zeros(np.shape(data_fe))
for mod in range(info_fe["Nmods"]):
    for lamb in range(info_fe["Nlambda"]):
        corrected_fe[0, lamb, mod, mask] = (data_fe[0, lamb, mod, mask]) / ff2[0, lamb, mod, mask]
        corrected_fe[1, lamb, mod, mask] = (data_fe[1, lamb, mod, mask]) / ff2[1, lamb, mod, mask]

hdu = fits.PrimaryHDU(corrected_mg)
hdu.writeto("obs_mode_1_ms.fits", overwrite = True)

hdu = fits.PrimaryHDU(corrected_fe)
hdu.writeto("obs_mode_2_ms.fits", overwrite = True)
"""
corrected_mg = fits.getdata("obs_mode_1_ms.fits")
corrected_fe = fits.getdata("obs_mode_2_ms.fits")


fig, axs = plt.subplots(2, 4, figsize = (20, 10))

axs[0, 0].set_title("OM 1 - core - cam1")
im = axs[0, 0].imshow(corrected_mg[0, 4, 0, 300:-300, 300:-300], cmap = "gray")
divider = make_axes_locatable(axs[0, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[0, 1].set_title("OM 1 - core - cam2")
im = axs[0, 1].imshow(corrected_mg[1, 4, 0, 300:-300, 300:-300], cmap = "gray")
divider = make_axes_locatable(axs[0, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)


axs[0, 2].set_title("OM 1 - cont - cam1")
im = axs[0, 2].imshow(corrected_mg[0, -1, 0, 300:-300, 300:-300], cmap = "gray")
divider = make_axes_locatable(axs[0, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[0, 3].set_title("OM 1 - cont - cam2")
im = axs[0, 3].imshow(corrected_mg[1, -1, 0, 300:-300, 300:-300], cmap = "gray")
divider = make_axes_locatable(axs[0, 3])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)


axs[1, 0].set_title("OM 2.02 - core - cam1")
im = axs[1, 0].imshow(corrected_fe[0, 3, 0, 300:-300, 300:-300], cmap = "gray")
divider = make_axes_locatable(axs[1, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[1, 1].set_title("OM 2.02 - core - cam2")
im = axs[1, 1].imshow(corrected_fe[1, 3, 0, 300:-300, 300:-300], cmap = "gray")
divider = make_axes_locatable(axs[1, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[1, 2].set_title("OM 1 - cont - cam1")
im = axs[1, 2].imshow(corrected_fe[0, -1, 0, 300:-300, 300:-300], cmap = "gray")
divider = make_axes_locatable(axs[1, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

axs[1, 3].set_title("OM 1 - cont - cam2")
im = axs[1, 3].imshow(corrected_fe[1, -1, 0, 300:-300, 300:-300], cmap = "gray")
divider = make_axes_locatable(axs[1, 3])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)


for i in range(2):
    for j in range(4):
        axs[i, j].set_xticks([])
        axs[i, j].set_yticks([])

plt.tight_layout()
plt.savefig(f"ms_presentation/continuum_corrected.png", bbox_inches = "tight")

demod, dual = demodulate(corrected_fe[:, -1], 2016, 2016, cf.om_config["2.02"]["Nmods"], cf.om_config["2.02"]["Nlambda"], filt = "525.02", mode = 'standard_single_wavelength') 

fig, axs = plt.subplots(3, 4, figsize = (14, 9))

norm = np.median(demod[0, 0, 300:-300, 300:-300])

axs[0, 0].set_title("I")
axs[0, 1].set_title("Q")
axs[0, 2].set_title("U")
axs[0, 3].set_title("V")

axs[0, 0].set_ylabel("Cam 1")
axs[1, 0].set_ylabel("Cam 2")
axs[2, 0].set_ylabel("Dual beam")


im = axs[0, 0].imshow(demod[0, 0, 200:1750, 200:1750] / norm, cmap = "gray", vmin = 0.8, vmax = 1.4)
divider = make_axes_locatable(axs[0, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[0, 1].imshow(demod[0, 1, 200:1750, 200:1750] / norm, cmap = "gray", vmin = -0.03, vmax = 0.03)
divider = make_axes_locatable(axs[0, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[0, 2].imshow(demod[0, 2, 200:1750, 200:1750] / norm, cmap = "gray", vmin = -0.03, vmax = 0.03)
divider = make_axes_locatable(axs[0, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[0, 3].imshow(demod[0, 3, 200:1750, 200:1750] / norm, cmap = "gray", vmin = -0.03, vmax = 0.03)
divider = make_axes_locatable(axs[0, 3])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

im = axs[1, 0].imshow(demod[1, 0, 200:1750, 200:1750] / norm, cmap = "gray", vmin = 0.8, vmax = 1.4)
divider = make_axes_locatable(axs[1, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[1, 1].imshow(demod[1, 1, 200:1750, 200:1750] / norm, cmap = "gray", vmin = -0.03, vmax = 0.03)
divider = make_axes_locatable(axs[1, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[1, 2].imshow(demod[1, 2, 200:1750, 200:1750] / norm, cmap = "gray", vmin = -0.03, vmax = 0.03)
divider = make_axes_locatable(axs[1, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[1, 3].imshow(demod[1, 3, 200:1750, 200:1750] / norm, cmap = "gray", vmin = -0.03, vmax = 0.03)
divider = make_axes_locatable(axs[1, 3])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

im = axs[2, 0].imshow(dual[0, 200:1750, 200:1750] / norm, cmap = "gray", vmin = 0.8, vmax = 1.4)
divider = make_axes_locatable(axs[2, 0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[2, 1].imshow(dual[1, 200:1750, 200:1750] / norm, cmap = "gray", vmin = -0.03, vmax = 0.03)
divider = make_axes_locatable(axs[2, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[2, 2].imshow(dual[2, 200:1750, 200:1750] / norm, cmap = "gray", vmin = -0.03, vmax = 0.03)
divider = make_axes_locatable(axs[2, 2])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
im = axs[2, 3].imshow(dual[3, 200:1750, 200:1750] / norm, cmap = "gray", vmin = -0.03, vmax = 0.03)
divider = make_axes_locatable(axs[2, 3])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)


for i in range(3):
    for j in range(4):
        axs[i, j].set_xticks([])
        axs[i, j].set_yticks([])

plt.tight_layout()
plt.savefig(f"ms_presentation/iron_demodulated.png", bbox_inches = "tight")


fig, axs = plt.subplots(figsize = (5.5, 5))
im = axs.imshow(dual[1, 200:1750, 200:1750] / norm, cmap = "gray")
divider = make_axes_locatable(axs)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.tight_layout()
plt.savefig(f"ms_presentation/iron_q_alone.png", bbox_inches = "tight")

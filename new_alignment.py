# ============================= IMPORTS ===================================== #
import os
import sys
import glob
import numpy as np 
import matplotlib.pyplot as plt

# Own Libs
sys.path.append("/home/users/dss/psanta/TuMags_Reduction_Pipeline")
import config as cf
from utils import read_Tumag
from field_stop_finder import compute_alignment, apply_fieldstop_and_align_array
from master_dark import compute_master_darks
from master_flatfield import compute_master_flat_field
import image_handler as ih
from demodulation import demodulate, demodulate_v2
import alignment as al

# ============================= CONFIG ===================================== #

dc_paths = ih.get_images_paths("D10-5620-5719")
ff_paths = ih.get_images_paths("D10-4340-5619")
Obs_paths = ih.get_images_paths("D10-304-1000")

# ======================= Darks and flats ================================== #

dc = compute_master_darks(dc_paths[:-4], verbose = True)
ff_data, ff_info = compute_master_flat_field(ff_paths,dc,  verbose = True)

# ======================= Selecting Observation Mode ======================= #

Ocs = ih.separate_ocs(Obs_paths, verbose = True, flat_fieldmode= False)
om = ih.nominal_observation("2.02", Ocs[84]["ims"], dc)
om_data = om.get_data()
om_info = om.get_info()

# ====================== Cropping and flat fielding ========================= #

# The 300 pixels remove the fieldstop
ff_cropped = ff_data[:, :, :, 300:-300, 300:-300]
om_data_cropped = om_data[:, :, :, 300:-300, 300:-300]

del ff_data, om_data # Deleting full arrays to save memory

# Flat-fielding
om_corr = np.zeros(np.shape(om_data_cropped))
for mod in range(om_info["Nmods"]):
    for lamb in range(om_info["Nlambda"]):
        om_corr[0, lamb, mod] = (om_data_cropped[0, lamb, mod]) / ff_cropped[0, lamb, mod]
        om_corr[1, lamb, mod] = (om_data_cropped[1, lamb, mod]) / ff_cropped[1, lamb, mod]


# =========================== Alignment and demodulation ====================== #

realigned = al.align_obsmode(om_corr, acc = 0.01, verbose = True)

dual, demod = demodulate_v2(realigned, om_info["Nmods"], om_info["Nlambda"], 
                            filt = om_info["line"], verbose = True )
# ---------------------------- DESCRIPTION --------------------------------------- #

"""
This module is an example for processing observing modes blocks. 

Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in Libs
import os
import sys
import numpy as np
from astropy.io import fits
import logging
import pickle
import pandas as pd

sys.path.append("/home/users/dss/psanta/TuMags_Reduction_Pipeline/")
#import main_processing as mp

# ------------------------------  CONFIG  ------------------------------------------ # 

Sunrise_id = "02_EFEM"
obs_block_idx = 0
observation_id = f"{Sunrise_id}_TM_{obs_block_idx}"
output_folder =  "/work/obs/TuMAG_data/Fits"
timeline_reduction_csv = "02_EMEF_data_processing_status.csv"

# ---------------------------------------------------------------------------- #
"""
dc = fits.getdata("/work/obs/TuMAG_data/calibration/darks_fits/dc_D10-6844-6943.fits")

ff_om1_1 = fits.getdata("/home/users/dss/psanta/Fits_creation/EMF_1/ff_om1_D10-36653-37236_D11-0-1655_avg_norm_pref_removed.fits")
ff_om1_2 = fits.getdata("/home/users/dss/psanta/Fits_creation/EMF_1/ff_om1_D10-7398-9637_avg_norm_pref_removed.fits")

ff_om1 = (ff_om1_1 + ff_om1_2) / 2
del ff_om1_1, ff_om1_2

# -------------------------------

ff_om202_1 = fits.getdata("/home/users/dss/psanta/Fits_creation/EMF_1/ff_om202_D10-9638-11429_avg_norm_pref_removed.fits")
ff_om202_2 = fits.getdata("/home/users/dss/psanta/Fits_creation/EMF_1/ff_om202_D11-1656-3447_avg_norm_pref_removed.fits")

ff_om202 = (ff_om202_1 + ff_om202_2) / 2
del ff_om202_1, ff_om202_2
"""
# ----------------------------------------------------------------------------

Compute_reduction_levels  = {
    "0.5" : True,
    "0.6" : True,
    "0.8" : True,
    "0.9" : True,
    "1.0" : True,
    "1.1" : True
}

Pickles = {
    "AR1" :  "ocs_ar1.pickle",
    "AR2" : "ocs_ar2.pickle"
}

# Configure logging
logging.basicConfig(
    filename = f"{obs_block_idx}_computation.log",  # Name of the log file
    level=logging.INFO,       # Log level (INFO, DEBUG, WARNING, ERROR, CRITICAL)
    format="%(asctime)s - %(levelname)s - %(message)s"  # Log format
)

def search_rows(csv_filename, conditions):
    df = pd.read_csv(csv_filename)

    # Apply multiple conditions
    query = " & ".join([f"{col} == @value" for col, value in conditions.items()])
    matching_rows = df.query(query, local_dict=conditions)

    return matching_rows



rows = search_rows(timeline_reduction_csv, {"status" : "Complete", "Obs_block" : "AR2"})
print(rows)

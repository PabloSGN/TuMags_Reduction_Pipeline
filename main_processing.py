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
import config
import image_handler as ih
from master_flatfield import correct_observation
from fits_files_handling import generate_fits
import phase_diversity as phase
import alignment as al
import demodulation as dem
import xtalk_jaeggli as jgg

# ------------------------------  CONFIG  ------------------------------------------ # 

Sunrise_id = "04_QSUN"
obs_block_idx = 0
observation_id = f"{Sunrise_id}_TM_{obs_block_idx}"
output_folder =  "/work/obs/TuMAG_data/Fits"

# Field-Stop
fs = [[155, 1799], 
      [171, 1815]]

lines = {"525.02" : "Fe",
         "525.06" : "Fe",
         "517" : "Mg"}

# ------------------------------  CODE  ------------------------------------------ # 

def update_csv_by_property(csv_filename, search_column, search_value, update_data):
    """
    Updates multiple column values for a row where a given column matches a specific value.

    Parameters:
        csv_filename (str): Path to the CSV file.
        search_column (str): Column to search for a specific value.
        search_value: The value to find in the search_column.
        update_data (dict): Dictionary with column names as keys and new values.
    """
    if not os.path.exists(csv_filename):
        raise FileNotFoundError(f"CSV file '{csv_filename}' not found.")

    df = pd.read_csv(csv_filename)

    if search_column not in df.columns:
        raise KeyError(f"Column '{search_column}' not found in CSV.")

    if search_value not in df[search_column].values:
        raise ValueError(f"Value '{search_value}' not found in column '{search_column}'.")

    for col in update_data.keys():
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in CSV.")

    df.loc[df[search_column] == search_value, update_data.keys()] = update_data.values()
    df.to_csv(csv_filename, index=False)


def compute_level_05(images_paths, observing_mode, dark_current, flat_field, field_stop):

    print("Processing reduction level 0.5...")

    om = ih.nominal_observation(observing_mode, images_paths, dark_current) 
    om_data = om.get_data()
    om_info = om.get_info()

    # Get the time for the filename
    t_start = om_info["Images_headers"][f"wv_0"]["M0"]["Date"].strftime('%d%m%YT%H%M%S')
    obsname_base = f"{observation_id}_{lines[om_info['line']]}{observing_mode}_{om_info['Nlambda']}_{t_start}"

    # Remove field-stop
    om_crop = om_data[:, :, :, field_stop[0][0]:field_stop[0][1], field_stop[1][0]:field_stop[1][1]]
    ff_crop = flat_field[:, :, :, field_stop[0][0]:field_stop[0][1], field_stop[1][0]:field_stop[1][1]]
    del om_data # Save memory

    om_corr = correct_observation(data = om_crop, ff = ff_crop)

    fits_file =  generate_fits(data = om_corr, fits_folder = output_folder, filename = obsname_base, 
                               level = 0.5, pipeline_version=config.Pipeline_version, om_info=om_info)
    
    return fits_file


def compute_level_06(flat_fielded_data, header, obsname, obsmode = None, zkes = None):

    print(f"Processing level 0.6 ...")

    if zkes is None and obsmode is not None:
        line = config.om_config[obsmode]["line"]
        zkes_id = f"{Sunrise_id}_{line}{obsmode}"
        reconstructed, zkes = phase.apply_reconstruction(data = flat_fielded_data, ID = zkes_id)
    elif zkes is not None:
        reconstructed, zkes = phase.apply_reconstruction(data = flat_fielded_data, zkes = zkes)
    else:
        raise Exception("If no zernikes are passed please provide obsmode")
    
    fits_file = generate_fits(data = reconstructed, fits_folder=output_folder, filename=obsname,
                              level = 0.6, pipeline_version=config.Pipeline_version, header = header)

    return fits_file

def compute_level_08_and_09(data, header, obsname, reductionlevel, acc = 0.01):

    aligned, shifts = al.align_obsmode(data, acc = acc, returnshifts = True)

    fits_file = generate_fits(data = aligned, fits_folder=output_folder, filename=obsname,
                           level = reductionlevel, pipeline_version=config.Pipeline_version,
                           header = header, shifts = shifts)

    return fits_file

def compute_level_1_and_higher(aligned_data, header, obsname, reductionlevel, shifts):

    line = header["FW2"]
    demod = dem.demodulate(aligned_data, filt = line)
    xtalked, mmatrix = jgg.fit_mueller_matrix(demod)

    fits_file = generate_fits(data = xtalked, fits_folder=output_folder, filename=obsname,
                           level = reductionlevel, pipeline_version=config.Pipeline_version,
                           header = header, shifts = shifts, fitted_muller=mmatrix)

    return fits_file

## ========================================================================= ##

def process_observation(Compute_reduction_levels, csvfile, Tstart, picklefile = None, oc = None, dc = None, ff = None, zkes = None, fits_file_05 = None,
                        fits_file_06 = None, fits_file_08 = None, fits_file_09 = None):

    logging.info(f"Starting reduction process..")
    try:
        # Computing level 0.5
        if Compute_reduction_levels["L0.5"]:

            with open(picklefile, "rb") as file:
                OCS = pickle.load(file) 

            images_paths = OCS[oc]["ims"]
            ObsMode = OCS[oc]["OM"]

            fits_file_05 = compute_level_05(images_paths=images_paths, observing_mode=ObsMode, 
                                            dark_current=dc, flat_field=ff)
            
            update_csv_by_property(csvfile, "Tstart", Tstart, {f"L0.5" : 1, "Pipeline" : config.Pipeline_version})
            logging.info(f"Level 0.5 - succesfully processed and stored. -> {fits_file_05}")

        if Compute_reduction_levels["L0.6"]:

            if fits_file_05 is not None:
                om_corr = fits.getdata(fits_file_05)
                header_05 = fits.getheader(fits_file_05)
                obsname_base = "_".join(os.path.basename(fits_file_05).split("_")[:7])
                ObsMode = header_05["OBS_MODE"]
            else:
                raise Exception("Provide path for level 0.5 fits or enable its computation")

            fits_file_06 = compute_level_06(flat_fielded_data = om_corr, header = header_05, 
                                            obsname = obsname_base, obsmode = ObsMode, zkes = zkes)
            update_csv_by_property(csvfile, "Tstart", Tstart, {f"L0.6" : 1, "Pipeline" : config.Pipeline_version})
            logging.info(f"Level 0.6 - succesfully processed and stored. -> {fits_file_06}")

        if Compute_reduction_levels["L0.8"]:
            
            if fits_file_05 is not None:
                om_corr = fits.getdata(fits_file_05)
                header_05 = fits.getheader(fits_file_05)
                obsname_base = "_".join(os.path.basename(fits_file_05).split("_")[:7])
                ObsMode = header_05["OBS_MODE"]
            else:
                raise Exception("Provide path for level 0.5 fits or enable its computation")

            fits_file_08 = compute_level_08_and_09(data = om_corr, header=header_05, 
                                                obsname=obsname_base, reductionlevel=0.8)
            
            update_csv_by_property(csvfile, "Tstart", Tstart, {f"L0.8" : 1, "Pipeline" : config.Pipeline_version})    
            logging.info(f"Level 0.8 - succesfully processed and stored. -> {fits_file_08}")
            
        if Compute_reduction_levels["L0.9"]:
            
            if fits_file_06 is not None:
                reconstructed = fits.getdata(fits_file_06)
                header_06 = fits.getheader(fits_file_06)
                obsname_base = "_".join(os.path.basename(fits_file_06).split("_")[:7])
                ObsMode = header_06["OBS_MODE"]
            else:
                raise Exception("Provide path for level 0.6 fits or enable its computation")

            fits_file_09 = compute_level_08_and_09(data = reconstructed, header=header_06, 
                                                obsname=obsname_base, reductionlevel=0.9)
            update_csv_by_property(csvfile, "Tstart", Tstart, {f"L0.9" : 1, "Pipeline" : config.Pipeline_version})
            logging.info(f"Level 0.9 - succesfully processed and stored. -> {fits_file_09}")
            
        if Compute_reduction_levels["L1.0"]:
            
            if fits_file_08 is not None:
                aligned = fits.getdata(fits_file_08)
                shifts = fits.getdata(fits_file_08, ext=1)
                header_08 = fits.getheader(fits_file_08)
                obsname_base = "_".join(os.path.basename(fits_file_08).split("_")[:7])
                ObsMode = header_08["OBS_MODE"]
            else:
                raise Exception("Provide path for level 0.8 fits or enable its computation")

            update_csv_by_property(csvfile, "Tstart", Tstart, {f"L1.0" : 1, "Pipeline" : config.Pipeline_version})
            fits_file_10 = compute_level_1_and_higher(aligned_data=aligned, header=header_08, 
                                                    obsname=obsname_base, reductionlevel=1.0, shifts=shifts)
            logging.info(f"Level 1.0 - succesfully processed and stored -> {fits_file_10}")
            
        if Compute_reduction_levels["L1.1"]:
            
            if fits_file_09 is not None:
                aligned = fits.getdata(fits_file_09)
                shifts = fits.getdata(fits_file_09, ext=1)
                header_09 = fits.getheader(fits_file_09)
                obsname_base = "_".join(os.path.basename(fits_file_09).split("_")[:7])
                ObsMode = header_09["OBS_MODE"]
            else:
                raise Exception("Provide path for level 0.8 fits or enable its computation")

            fits_file_11 = compute_level_1_and_higher(aligned_data=aligned, header=header_09, 
                                                    obsname=obsname_base, reductionlevel=1.1, shifts=shifts)
            
            update_csv_by_property(csvfile, "Tstart", Tstart, {f"L1.1" : 1, "Pipeline" : config.Pipeline_version})
            logging.info(f"Level 1.1 - succesfully processed and stored -> {fits_file_11}")

    except Exception as e:
        logging.error(f"ERROR in processing OC -> {oc}. Error -> {e}")

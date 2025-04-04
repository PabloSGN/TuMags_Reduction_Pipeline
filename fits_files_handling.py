# ---------------------------- DESCRIPTION --------------------------------------- #

"""
This module provides functionality to handle and process fits files. 

Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in Libs
import os
import sys
import numpy as np
from astropy.io import fits

# ------------------------------  CONFIG  ------------------------------------------ # 

spectral_cal = {
'517' : {
         'Pend'     : 0.00030907042253499933,
         'Ord'      : 5173.432608450703},
'525.02': {
        'Pend'     : 0.0002957121398329138,
        'Ord'      : 5249.543594995222},
'525.06' : {
        'Pend'     : 0.000288733333332857,
        'Ord'      : 5251.371833333332}}

# ------------------------------  CODE  ------------------------------------------ # 

def volts_2_lambda(volts, config):
    return config['Pend'] * volts + config['Ord']

def header_generator(ominfo):

    """
    Function that generates the header for the fits file from the observing 
    modes information -> om.get_info()
    Inputs:
      - ominfo : Dictionary returned by the om.get_info() routine. 
    Returns:
      - Header for fits files. 
      - comms : Comments of header parameters
    """

    h = {
    "Obs_Mode" : ominfo["ObservationMode"],
    "nAcc" : ominfo["nAcc"],
    "Roix" : ominfo["Roix"],
    "Roiy" : ominfo["Roiy"],
    "Nlambda" : ominfo["Nlambda"],
    "T_exp" : ominfo["Images_headers"][f"wv_0"]["M0"]["t_exp"],
    "FW1" : ominfo["Images_headers"][f"wv_0"]["M0"]["FW1"],
    "FW2" : ominfo["Images_headers"][f"wv_0"]["M0"]["FW2"],
    "Nmods" : ominfo["Nmods"],
    "OC" : ominfo["Images_headers"][f"wv_0"]["M0"]["ObservationCounter"],
    "t_START" : ominfo["Images_headers"][f"wv_0"]["M0"]["Date"].strftime("%d/%m/%Y, %H:%M:%S"),
    "t_END" : ominfo["Images_headers"][f"wv_{ominfo['Nlambda'] - 1}"]["M3"]["Date"].strftime("%d/%m/%Y, %H:%M:%S")
    }

    comms = {
    "Obs_Mode" : "Observation Mode",
    "nAcc" : "Number of accumulations",
    "Roix" : "Size in x",
    "Roiy" : "Size in y",
    "Nlambda" : "Number of wavelengths",
    "T_exp" : "Exposure time [ms]",
    "FW1" : "Filter wheel 1 position",
    "FW2" : "Spectral line (FW2 pos.)",
    "Nmods" : "Number of modulations",
    "OC" : "Observation Cunter",
    "t_START" : "Time of first image",
    "t_END" : "Time of last image"
    }

    for ind, ll in enumerate(ominfo["lambda_array"]):
        #h[f"LBD_{ind}"] = round(ll, 3)
        h[f"L_{ind}"] = round(volts_2_lambda(ominfo["Images_headers"][f"wv_{ind}"]["M0"]["hvps_read_volts"],
                                         spectral_cal[ominfo["line"]]), 3)
        comms[f"L_{ind}"] = f"Wavelength P{ind} [A]"
        
        h[f"V_{ind}"] = round(ominfo["Images_headers"][f"wv_{ind}"]["M0"]["hvps_read_volts"], 3)
        comms[f"V_{ind}"] = f"Etalon volts P{ind} [V]"

    return h, comms


def generate_fits(data, fits_folder, filename, level, pipeline_version, om_info = None, 
                  header = None, zkes = None,  shifts = None, fitted_muller = None):

    """
    Function that saves a fits file and generates header for a given observing mode. 
    Inputs:
       - data : Obs mode (np.array)
       - om_info : Dictionary returned by the om.get_info() routine. 
       - fits_folder : Path to store the data on. 
       - filename : Basename filename (without reduction level info)
       - level : Level of reduction. 
       - pipeline_version : Pipeline version
       - Zkes : Zernike coefficients used. 
       - shifts : Shifts applied to images in alignment
       - fitted_muller : Matrix fitted in xtalk analysis. 
    """

    # Create the filename with the reduction level
    extended_filename = f"{filename}_LV_{level}_v{pipeline_version}.fits"

    # If adding info to a lower reduction level with header already created
    if header is not None:
        hdu = fits.PrimaryHDU(data.astype(np.float32), header = header)
    
    # Creating the header from om info
    elif om_info is not None:
        # Create hdu
        hdu = fits.PrimaryHDU(data.astype(np.float32))
        header = hdu.header # Extract header
        h, comments = header_generator(om_info) # create header from om_info

        # Store om_info in HDU's header
        header["ppl_ver"] = pipeline_version 
        header.comments["ppl_ver"] = "Pipeline version"

        header["level"] = level
        header.comments["level"] = "Level of reduction"

        for key in h:
            header[key] = h[key]
            header.comments[key] = comments[key]

    else:
        raise Exception("Please provide om_info or already created header to generate_fits")

    # Store Zernikes if passed. (For level 0.6)
    if zkes is not None:
        for zkind, zk in enumerate(zkes):
            header[f"Z_{zkind}"] = round(zk, 4)
            header.comments[f"Z_{zkind}"] = f"Zernike index {zkind}"

    hdul_flag = False
    # Store alignment shifts if passed. (for levels 0.8 and 0.9)
    if shifts is not None:
        hdul_flag = True
        shifts_hdu = fits.ImageHDU(shifts.astype(np.float32), name="Shifts")

        if fitted_muller is None:
            hdul = fits.HDUList([hdu, shifts_hdu])

    # Store fitted mueller matrix if passed. (for levels 1.0 and 1.1)
    if fitted_muller is not None:
        hdul_flag = True
        muller_hdu = fits.ImageHDU(fitted_muller.astype(np.float32), name = "Fitted_Muller_Matrix")
        hdul = fits.HDUList([hdu, shifts_hdu, muller_hdu])

    # Save the file. 
    if hdul_flag:
        hdul.writeto(os.path.join(fits_folder, extended_filename), overwrite=True)
    else:
        hdu.writeto(os.path.join(fits_folder, extended_filename), overwrite=True)
        
    return os.path.join(fits_folder, extended_filename)


    

    

    

   
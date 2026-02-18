# ---------------------------- DESCRIPTION --------------------------------------- #

"""
This module provides functionality to handle and process fits files. 

Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Standard library
import os
from datetime import datetime

# Third-party libraries
import numpy as np
from astropy.io import fits
from astropy.time import Time

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
def update_header(header, keyword, value, after=None, comment=None):
    if keyword in header:  # Check for existence
        header[keyword] = value
    else:
        header.set(keyword, value, comment, after=after)
    return header

def volts_2_lambda(volts, config):
    return config['Pend'] * volts + config['Ord']

def header_generator(data, ominfo,datatype="SCIENCE"):

    """
    Function that generates the header for the fits file from the observing 
    modes information -> om.get_info()
    Inputs:
      - ominfo : Dictionary returned by the om.get_info() routine. 
    Returns:
      - Header for fits files. 
      - comms : Comments of header parameters
    """

    shape = data.shape
    ndim = data.ndim
    ny, nx = shape[-2], shape[-1]

    obs_start = ominfo["Images_headers"]["wv_0"]["M0"]["Date"]
    obs_end = ominfo["Images_headers"][f"wv_{ominfo['Nlambda'] - 1}"]["M3"]["Date"]

    h = {
        "SIMPLE"   : True,
        "BITPIX"   : -32,  # 32-bit float
        "NAXIS"    : ndim,
        "NAXIS1"   : ominfo["Roix"],
        "NAXIS2"   : ominfo["Roiy"],
        "EXTEND"   : True,
        "BUNIT"    : "DN",
        "DATATYPE" : datatype,  # DARK / FLAT / SCIENCE
        # Instrument
        "INSTRUME" : "TuMag",
        "TELESCOP" : "Sunrise",
        "OBSERVER" : "IAA-CSIC",
        # Observation timing
        "DATE-OBS" : obs_start.isoformat(),
        "DATE-END" : obs_end.isoformat(),
        "T_START"  : obs_start.strftime("%d/%m/%Y, %H:%M:%S"),
        "T_END"    : obs_end.strftime("%d/%m/%Y, %H:%M:%S"),
        # Image metadata
        "EXPTIME"  : ominfo["Images_headers"]["wv_0"]["M0"]["t_exp"],  # in ms
        "FW1"      : ominfo["Images_headers"]["wv_0"]["M0"]["FW1"],
        "FW2"      : ominfo["Images_headers"]["wv_0"]["M0"]["FW2"],
        "OBSMODE"  : ominfo["ObservationMode"],
        "NACC"     : ominfo["nAcc"],
        "NWAVE"    : ominfo["Nlambda"],
        "NMODS"    : ominfo["Nmods"],
        "OBSCNT"   : ominfo["Images_headers"]["wv_0"]["M0"]["ObservationCounter"],
        # Wavelength info (example values)
        "WAVE0"    : ominfo["Images_headers"]["wv_0"]["M0"]["FW2"],       
        "WAVEUNIT": "nm",
        # WCS Solar coordinates (HMI-style)
        "CTYPE1"   : "HPLN-TAN",
        "CTYPE2"   : "HPLT-TAN",
        "CUNIT1"   : "arcsec",
        "CUNIT2"   : "arcsec",
        "CRPIX1"   : ominfo["Roix"] // 2,
        "CRPIX2"   : ominfo["Roiy"] // 2,
        "CRVAL1"   : 0.0,
        "CRVAL2"   : 0.0,
        "CDELT1"   : 0.5,         # arcsec/pixel (adjust as needed)
        "CDELT2"   : 0.5,
        "SOLAR_P"  : 0.0,
        "RSUN_OBS" : 976.0,
        "DSUN_OBS" : 1.496e+11,
        # --- Calibration / processing provenance ---
        "LINECORR" : 0,
        "LINEPARS" : 0,
        "LINEMETH" : "",

        "DARKCORR" : 0,
        "DARK_ID"  : "",

        "FLATCORR" : 0,
        "FLAT_ID"  : "",

        "FILTERED" : 0,
        "FILTMETH" : "",
        "FILT_ID"  : "",

        "REALIGN"  : 0,
        "ALIGN_ID" : "",
        "ALIGMETH" : "",

        "WFR_RECO" : 0,
        "WFR_ID"   : "",
        "WFR_METH" : "",

        # --- versioning ---
        "PIPE_VER" : "",
        "PROCDATE" : Time.now().isot,
        "HISTORY"  : " ",
    }

    if ndim >= 3:
        h["NAXIS3"] = shape[-3]  # e.g. NWAVE
    if ndim == 4:
        h["NAXIS4"] = shape[-4]  # e.g. NMODS

    comms = {
        "OBSMODE"   : "Observation Mode",
        "NACC"      : "Number of accumulations",
        "NWAVE"     : "Number of wavelengths",
        "EXPTIME"   : "Exposure time [ms]",
        "FW1"       : "Filter wheel 1 position",
        "FW2"       : "Spectral line (FW2 pos.)",
        "NMODS"     : "Number of modulation states",
        "OBSCNT"    : "Observation Counter",
        "T_START"   : "Start time (first frame)",
        "T_END"     : "End time (last frame)",
        # Calibration status
        "LINECORR"  : "Was linearity correction applied?",
        "LINEPARS"  : "Linearity coefficients",
        "LINEMETH"  : "Linearity correction method",
        "DARKCORR"  : "Was dark correction applied?",
        "DARK_ID"   : "Dark frame ID",
        "FLATCORR"  : "Was flatfielding applied?",
        "FLAT_ID"   : "Flat frame ID",
        "FILTERED"  : "Was filtering applied?",
        "FILT_ID"   : "Filter version or ID",
        "FILTMETH"  : "Filter method",
        "REALIGN"   : "Was alignment applied?",
        "ALIGN_ID"  : "Alignment session or method ID",
        "ALIGMETH"  : "Alignment method",
        "WFR_RECO"  : "Was wavefront reconstruction applied?",
        "WFR_ID"    : "WFR configuration ID",
        "WFR_METH"  : "Wavefront reconstruction method",
        # Version tracking
        "PIPE_VER"  : "Reduction pipeline version",
        "PROCDATE"  : "Date/time of calibration",
        "HISTORY"   : "Processing history summary",
        # Wavelength info
        "WAVE0"     : "Base wavelength [nm]",
        "WAVEUNIT"  : "Wavelength unit",
        # WCS Solar coordinates (HMI-style)
        "CTYPE1"    : "WCS axis 1 type (solar-x)",
        "CTYPE2"    : "WCS axis 2 type (solar-y)",
        "CUNIT1"    : "WCS axis 1 unit",
        "CUNIT2"    : "WCS axis 2 unit",
        "CRPIX1"    : "WCS reference pixel X",
        "CRPIX2"    : "WCS reference pixel Y",
        "CRVAL1"    : "WCS reference value X",
        "CRVAL2"    : "WCS reference value Y",
        "CDELT1"    : "WCS pixel scale X [arcsec/pix]",
        "CDELT2"    : "WCS pixel scale Y [arcsec/pix]",
        "SOLAR_P"   : "Solar P angle [deg]",
        "RSUN_OBS"  : "Observed solar radius [arcsec]",
        "DSUN_OBS"  : "Sun-observer distance [m]",
        "SIMPLE"    : "FITS: file conforms to standard",
        "BITPIX"    : "FITS: number of bits per data pixel",
        "NAXIS"     : "FITS: number of data axes",
        "NAXIS1"    : "FITS: axis 1 length",
        "NAXIS2"    : "FITS: axis 2 length",
        "NAXIS3"    : "FITS: axis 2 length",
        "NAXIS4"    : "FITS: axis 2 length",
        "EXTEND"    : "FITS: may contain extensions",
        "BUNIT"     : "Physical units of array values",
        "DATATYPE"  : "Data type: DARK / FLAT / SCIENCE",
        "INSTRUME"  : "Instrument name",
        "TELESCOP"  : "Telescope name",
        "OBSERVER"  : "Observer / institution",
        "DATE-OBS"  : "Observation start (ISO)",
        "DATE-END"  : "Observation end (ISO)",
    }


    for ind, ll in enumerate(ominfo["lambda_array"]):
        #h[f"LBD_{ind}"] = round(ll, 3)
        h[f"L_{ind}"] = round(volts_2_lambda(ominfo["Images_headers"][f"wv_{ind}"]["M0"]["hvps_read_volts"],
                                         spectral_cal[ominfo["line"]]), 3)
        comms[f"L_{ind}"] = f"Wavelength P{ind} [A]"
        
        h[f"V_{ind}"] = round(ominfo["Images_headers"][f"wv_{ind}"]["M0"]["hvps_read_volts"], 3)
        comms[f"V_{ind}"] = f"Etalon volts P{ind} [V]"

    for wl in ominfo["Images_headers"]:
        for pn in ominfo["Images_headers"][wl]:
            file_id = ominfo["Images_headers"][wl][pn]['image_name']
            # print(file_id)
            dt_str = file_id.split('_0_')[0]  # '2024_07_12_09_27_29_344'
            dt = datetime.strptime(dt_str, '%Y_%m_%d_%H_%M_%S_%f')            
            h[f"{wl}_{pn}"] = dt.strftime("%d%m%YT%H%M%S.%f")
            comms[f"{wl}_{pn}"] = f"Exact wave adquisition"

    for wl in ominfo["Images_headers"]:
        for pn in ominfo["Images_headers"][wl]:
            lcvr1_volts = ominfo["Images_headers"][wl][pn]['lcvr1_volts']
            lcvr2_volts = ominfo["Images_headers"][wl][pn]['lcvr2_volts']
            h[f"L1{wl}{pn}"] = lcvr1_volts
            h[f"L2{wl}{pn}"] = lcvr2_volts
            comms[f"L1{wl}{pn}"] = f"Read lcvr1 in volts"
            comms[f"L2{wl}{pn}"] = f"Read lcvr1 in volts"

    return h, comms


def generate_fits(data, fits_folder, extended_filename, level, pipeline_version, om_info = None, 
                  header = None, zkes = None,  shifts = None, fitted_muller = None,datatype=None,
                  FLAT_ID = None, DARK_ID = None,
                  ):

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
    # extended_filename = f"{filename}_LV_{level}_v{pipeline_version}.fits"

    # If adding info to a lower reduction level with header already created
    if header is not None:
        hdu = fits.PrimaryHDU(data.astype(np.float32), header = header)
    
    # Creating the header from om info
    elif om_info is not None:
        # Create hdu
        hdu = fits.PrimaryHDU(data.astype(np.float32))
        header = hdu.header # Extract header
        h, comments = header_generator(data, om_info,datatype=datatype) # create header from om_info

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

    # Store fitted mueller matrix if passed. (for levels 1.0 and 1.1)
    if FLAT_ID is not None:
        header["FLATCORR"] = 1
        header["FLAT_ID"] = FLAT_ID
        
    if DARK_ID is not None:
        header["DARKCORR"] = 1
        header["DARK_ID"] = DARK_ID
        


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


    # add additional keyword elements: 
    # Save the file. 
    if hdul_flag:
        hdul.writeto(os.path.join(fits_folder, extended_filename), overwrite=True)
    else:
        hdu.writeto(os.path.join(fits_folder, extended_filename), overwrite=True)
        
    return os.path.join(fits_folder, extended_filename)


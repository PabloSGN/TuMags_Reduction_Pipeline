import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import pandas as pd


def get_z_elevation(elevation,nz=6,filter='517.3'):
    
    #Import numpy array with the fit parameters (slope and intercept) for each Zernike coefficient
    #  filter='517.3' #'525.02', '525.06' or '517.3' #Filter to be used
    
    MODULE_DIR = Path(__file__).resolve().parent 

    filepath = MODULE_DIR / "CALDATA" / f"fit_zernike_params_{filter}.npy"

    fit_params=np.load(filepath)

    Nzer=fit_params.shape[0] #Number of Zernike coefficients

    #Use the fit parameters to model the Zernike coefficients
    # elevation=25 #Elevation in degrees


    zernikes_model=np.zeros((Nzer+3))
    # for k in range(3,Nzer+3):
    for k in range(3,nz+3):
        zernikes_model[k]=fit_params[k-3,0]*elevation+fit_params[k-3,1]

    return zernikes_model

def get_z(filter='517.3'):
    
    #Import numpy array with the fit parameters (slope and intercept) for each Zernike coefficient
    #  filter='517.3' #'525.02', '525.06' or '517.3' #Filter to be used
    MODULE_DIR = Path(__file__).resolve().parent

    filepath = MODULE_DIR / "CALDATA" / f"fit_zernike_params_{filter}.npy"

    print(filepath)  # para depurar

    return np.load(filepath)

# def get_elevation():

#     MODULE_DIR = Path(__file__).resolve().parent

#     file = MODULE_DIR / "CALDATA" / f'Sun_elevation_pointing_data.csv'

#     data = pd.read_csv(
#         file,
#         delimiter=";",
#         usecols=["Time", "Sun elevation [deg]"]
#     )

#     data["Time"] = pd.to_datetime(data["Time"])

#     return data

def get_elevation(header):

    MODULE_DIR = Path(__file__).resolve().parent

    file = MODULE_DIR / "CALDATA" / f'Sun_elevation_pointing_data.csv'

    elevation_data = pd.read_csv(
        file,
        delimiter=";",
        usecols=["Time", "Sun elevation [deg]"]
    )

    elevation_data["Time"] = pd.to_datetime(
        elevation_data["Time"],
        format="mixed"
    )
    obs_time = pd.to_datetime(header["DATE-OBS"])

    idx = (elevation_data["Time"] - obs_time).abs().idxmin()

    elevation = elevation_data.loc[idx, "Sun elevation [deg]"]

    return float(elevation)
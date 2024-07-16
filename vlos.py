# ---------------------------- DESCRIPTION --------------------------------------- #

# ------------------------------ IMPORTS ----------------------------------------- #
import numpy as np
import time

from utils import read_Tumag

from scipy.integrate import simps

# ------------------------------ CONFIG ------------------------------------------ #

# ------------------------------  CODE  ------------------------------------------ # 


def cog_profile(input_data, wave_axis, wave, cpos = -1):


    input_data = np.transpose(input_data, axes = (1, 2, 0))

    Ic = np.max(input_data, axis = -1)

    tc = []
    for i in range(len(wave_axis)):
        tc.append(( Ic - input_data[:, :, i] ))

    tc = np.transpose(np.array(tc), axes = (1, 2, 0))

    lcog = simps(wave_axis * tc, x=wave_axis, axis = -1) / simps(tc, x=wave_axis, axis = -1)

    vlos = - (wave - lcog ) * 2.99792458e+5 / wave

    return vlos

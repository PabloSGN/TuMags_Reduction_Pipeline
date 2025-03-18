# ---------------------------- DESCRIPTION --------------------------------------- #

"""
Demodulation functions
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in libs
import numpy as np

# Own-libs
import alignment as al
# ------------------------------ CONFIG ------------------------------------------ #

# Mean - Matrix Demodulation 
mod_matrices = { # Calculadas por Antonio C -> 18 Abril Kiruna 2024 
    "517": {0 : np.array([[0.951,  -0.612,	0.474	,0.459],
                          [0.955,	-0.331,	-0.758,	-0.382],
                          [1.058,	 0.456,	0.562   ,-0.712],
                          [1.036,	 0.747,	-0.26	,0.6]]),

            1 : np.array([[1.054,	0.763,  -0.394, -0.524],
                          [1.036,	0.497,	0.793,	 0.306],
                          [0.953,  -0.282,  -0.475,	 0.683],
                          [0.958,  -0.585,	0.32,   -0.613]])},

    "525.02" : {0 : np.array([[0.954,	-0.694,	 0.406,  0.414],
                              [0.969,	-0.39,	-0.803,	-0.368],
                              [1.042,  0.418,	 0.495,	-0.705],
                              [1.035,  0.71,	-0.266,	 0.612 ]]),
                                
                1 : np.array([[1.059,	0.771,  -0.449,	-0.433],
                              [1.024,	0.449,   0.723,	 0.335],
                              [0.965,	-0.344,	-0.543,	 0.65],
                              [0.953,	-0.606,	 0.191,	-0.641]])},  

    "525.06" : {0 : np.array([[0.951,	-0.687,	 0.403,	 0.424],
                              [0.962,	-0.373,	-0.8,	-0.339],
                              [1.048,	 0.415,	 0.5,	-0.728],
                              [1.038,	 0.736,	-0.236,	 0.601]]), 

                1 : np.array([[1.06	 ,   0.777,	-0.403,	-0.463],
                              [1.032,    0.471,	 0.754,	 0.29],
                              [0.96	 ,  -0.306,	-0.497,	 0.681],
                              [0.948 ,  -0.62,	 0.205,	-0.619]])},            
 }


mod_matrices_david = { 
    "517": {0 : np.array([[0.9655, -0.4865,  0.6307, 0.4986],
                          [0.9476, -0.5615, -0.6319, -0.3653],
                          [1.0471,  0.5569,  0.4372, -0.7102],
                          [1.0398,  0.6294, -0.4595, 0.6237]]),

            1 : np.array([[1.0505, 0.5992 , -0.6032 ,-0.5262],
                          [1.0372, 0.6798 , 0.6194  ,0.3431],
                          [0.9663, -0.4213, -0.4031 , 0.7268],
                          [0.9459, -0.5206, 0.4653  ,-0.5974]])},

    "525.02" : {0 : np.array([[0.9603, -0.5244,  0.6005, 0.4470],
                              [0.9525, -0.5592, -0.6413, -0.3498],
                              [1.0435, 0.5823 ,  0.3865, -0.7071],
                              [1.0437, 0.6645 , -0.4272, 0.6306 ]]),
                                
                1 : np.array([[1.0598, 0.6811 , -0.6278, -0.4302],
                              [1.0408, 0.6922 , 0.5852 ,  0.3614],
                              [0.9610, -0.4278, -0.4182,  0.6859],
                              [0.9384, -0.4984, 0.3609 , -0.6310]])},  

    "525.06" : {0 : np.array([[0.9557, -0.5389,  0.5895, 0.4505],
                              [0.9486, -0.5501, -0.6593, -0.3247],
                              [1.0485, 0.5693 , 0.3855 , -0.7246],
                              [1.0473, 0.6737 , -0.4148, 0.6236]]), 

                1 : np.array([[1.0635, 0.6843 , -0.5910 ,-0.4638],
                              [1.0466, 0.7028 , 0.6070  ,0.3110],
                              [0.9570, -0.3960, -0.3864 ,0.7157],
                              [0.9329, -0.5181, 0.3745  ,-0.6051]])},            
 }


# Compute demodulation matrixes by inverting
demod_matrices = {}
for filt in mod_matrices:
    demod_matrices[filt] = {}
    for cam in mod_matrices[filt]:
        demod_matrices[filt][cam] = np.linalg.inv(mod_matrices[filt][cam])

demod_matrices_david = {}
for filt in mod_matrices_david:
    demod_matrices_david[filt] = {}
    for cam in mod_matrices_david[filt]:
        demod_matrices_david[filt][cam] = np.linalg.inv(mod_matrices_david[filt][cam])

# ------------------------------  CODE  ------------------------------------------ # 

def demodulate(data, filt, dmod_matrices = demod_matrices_david, onelambda = False, BothCams = False):
    """
    Function to perform the demodulation of the observation mode. 
    Inputs: 
        - data (np.array) : Array contaning the obs mode. (Ncams x Nlambda x Nmods x Nx x Ny)
        - filt (str) : Filter to demodulate (517, 525.02 or 525.06)
        - demod_matrices (np.array : default : demod_matrices_David) : Demodulations matrix to use
       - onelambda (Boolean, default : False): Set to true if only one lambda is used (array of shape Ncam x Nmod x Nx x Ny) 
    Outputs:
        - dual_beamed (np.array) : Demodulated data with cameras combined (Nlambda x Nmods x Nx x Ny).
        - demodulated (np.array) : Demodulated data with cameras not yet combined (Ncams x Nlambda x Nmods x Nx x Ny). 
    """

    if onelambda:
        data = data[:, np.newaxis] # To allow for only one lamdba.

    shape = np.shape(data)
    nlambda = shape[1]
    nmods = shape[2]
    
    # All wavelengths
    size = np.shape(data)[-1]
    demod = np.zeros(np.shape(data))
    dual_beam = np.zeros((nlambda, nmods, size, size))

    # Each wavelength independently
    for wl in range(nlambda):
        
        dm_cam1 = np.matmul(dmod_matrices[filt][0], np.reshape(data[0, wl, :], (4, size * size)))
        dm_cam2 = np.matmul(dmod_matrices[filt][1], np.reshape(data[1, wl, :], (4, size * size)))

        demod[0, wl, :] = np.reshape(dm_cam1, (4, size, size))
        demod[1, wl, :] = np.reshape(dm_cam2, (4, size, size))
    
        dual_beam[wl] = (demod[0, wl] + demod[1, wl]) / 2
    
    if BothCams:
        if onelambda:
            return dual_beam[0], demod[:, 0]
        else:
            return dual_beam, demod
    else:
        if onelambda:
            return dual_beam[0]
        else:
            return dual_beam


def demodulate_quadrants(data, nlambda, nmods, filt, nquads = 16, Np_quad = 354):
    """
    Function to perform the demodulation of the observation mode separated in quadrants. 
    Inputs: 
        - data (np.array) : Array contaning the obs mode. (Ncams x Nlambda x Nmods x Nx x Ny)
        - nlambda (int): Number of wavelengths
        - nmods (int) : Number of modulations
        - filt (str) : Filter to demodulate (517, 525.02 or 525.06)
        - nquads (int, default : 16) : Number of quadrants
        - Np_quad (int, defaulr : 354) : Pixel soize of quadrant
    Outputs:
        - dual_beamed (np.array) : Demodulated data with cameras combined (Nlambda x Nmods x Nx x Ny).
        - demodulated (np.array) : Demodulated data with cameras not yet combined (Ncams x Nlambda x Nmods x Nx x Ny). 
    """

    demod = np.zeros((2, nlambda, nmods, nquads, Np_quad, Np_quad))
    dual = np.zeros((nlambda, nmods, nquads, Np_quad, Np_quad))

    for quad in range(nquads):

        du, dem = demodulate(data[:, :, :, quad], nmods, nlambda, filt)
        demod[:, :, :, quad] = dem
        dual[:, :, quad] = du

    return dual, demod

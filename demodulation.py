# ---------------------------- DESCRIPTION --------------------------------------- #

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in libs
import numpy as np

# Own-libs

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

def demodulate(data, sizex, sizey, nmods, nlambdas, filt, mode = 'standard', dmod_matrices = demod_matrices_david):
   
    if mode == 'standard':

        demod = np.zeros(np.shape(data))
        dual_beam = np.zeros((nlambdas, nmods, sizex, sizey))

        for wl in range(nlambdas):

            dm_cam1 = np.matmul(dmod_matrices[filt][0], np.reshape(data[0, wl, :], (4, sizex * sizey)))
            dm_cam2 = np.matmul(dmod_matrices[filt][1], np.reshape(data[1, wl, :], (4, sizex * sizey)))

            demod[0, wl, :] = np.reshape(dm_cam1, (4, sizex, sizey))
            demod[1, wl, :] = np.reshape(dm_cam2, (4, sizex, sizey))

        dual_beam = demod[0] * 0.5 + demod[1] * 0.5

    if mode == 'standard_single_wavelength':

        demod = np.zeros(np.shape(data))
        dual_beam = np.zeros((nmods, sizex, sizey))

        dm_cam1 = np.matmul(dmod_matrices[filt][0], np.reshape(data[0], (4, sizex * sizey)))
        dm_cam2 = np.matmul(dmod_matrices[filt][1], np.reshape(data[1], (4, sizex * sizey)))

        demod[0] = np.reshape(dm_cam1, (4, sizex, sizey))
        demod[1] = np.reshape(dm_cam2, (4, sizex, sizey))

        dual_beam = demod[0] * 0.5 + demod[1] * 0.5

    """
    elif mode == 'pixel':

        Npix = 1656
        data = np.zeros((2, 4, Npix, Npix))

        for ind, im in enumerate(sorted_files[filt_index, wvl_index, :, 0]):
            I, _ = read_Tumag(im)
            
            data[0, ind] = I[dmod_px2px_c1["RECT"][0][0]:dmod_px2px_c1["RECT"][0][0] + dmod_px2px_c1["RECT"][0][2] + 1,
                             dmod_px2px_c1["RECT"][0][1]:dmod_px2px_c1["RECT"][0][1] + dmod_px2px_c1["RECT"][0][3] + 1]

        for ind, im in enumerate(sorted_files[filt_index, wvl_index, :, 1]):
            I, _ = read_Tumag(im)
            data[1, ind] = I[dmod_px2px_c2["RECT"][0][0]:dmod_px2px_c2["RECT"][0][0] + dmod_px2px_c2["RECT"][0][2] + 1,
                             dmod_px2px_c2["RECT"][0][1]:dmod_px2px_c2["RECT"][0][1] + dmod_px2px_c2["RECT"][0][3] + 1]

        data[0] /= np.max(data[0])
        data[1] /= np.max(data[1])

        demod = np.zeros(np.shape(data))

        #d1 = np.transpose(dmod_px2px_c1["O"], axes = (2, 3, 0, 1))
        #d2 = np.transpose(dmod_px2px_c2["O"], axes = (2, 3, 0, 1))

        demod[0] = np.einsum('ijlk,kij->lij', np.transpose(dmod_px2px_c1["O"], axes =(3, 2)), data[0])
        demod[1] = np.einsum('ijlk,kij->lij', dmod_px2px_c2["O"], data[1]) 
        
        dual_beam = demod[0] * 0.5 + np.flip(demod[1], axis = 2) * 0.5

    elif mode == 'pixel_slow':

        print("Started slooooow dem...")
        Npix = 1656
        data = np.zeros((2, 4, Npix, Npix))

        for ind, im in enumerate(sorted_files[filt_index, wvl_index, :, 0]):
            I, _ = read_Tumag(im)
            
            data[0, ind] = I[dmod_px2px_c1["RECT"][0][0]:dmod_px2px_c1["RECT"][0][0] + dmod_px2px_c1["RECT"][0][2] + 1,
                             dmod_px2px_c1["RECT"][0][1]:dmod_px2px_c1["RECT"][0][1] + dmod_px2px_c1["RECT"][0][3] + 1]

        for ind, im in enumerate(sorted_files[filt_index, wvl_index, :, 1]):
            I, _ = read_Tumag(im)
            data[1, ind] = I[dmod_px2px_c2["RECT"][0][0]:dmod_px2px_c2["RECT"][0][0] + dmod_px2px_c2["RECT"][0][2] + 1,
                             dmod_px2px_c2["RECT"][0][1]:dmod_px2px_c2["RECT"][0][1] + dmod_px2px_c2["RECT"][0][3] + 1]

        data[0] /= np.max(data[0])
        data[1] /= np.max(data[1])

        demod = np.zeros(np.shape(data))

        for x in range(Npix):
            print(f"R:{x}/{Npix}")
            for y in range(Npix):
                
                dm_cam1 = np.matmul(dmod_px2px_c1["O"][y, x], np.reshape(data[0], (4, Npix * Npix)))
                dm_cam2 = np.matmul(dmod_px2px_c2["O"][y, x], np.reshape(data[1], (4, Npix * Npix)))

                demod[0] = np.reshape(dm_cam1, (4, Npix, Npix))
                demod[1] = np.reshape(dm_cam2, (4, Npix, Npix))

        dual_beam = demod[0] * 0.5 + np.flip(demod[1], axis = 2) * 0.5"""

    return demod, dual_beam

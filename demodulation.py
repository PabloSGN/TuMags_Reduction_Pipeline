# ---------------------------- DESCRIPTION --------------------------------------- #

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in libs
import numpy as np

# Own-libs

# ------------------------------ CONFIG ------------------------------------------ #

# Mean - Matrix Demodulation 
mod_matrices = { # Calculadas por Antonio C -> 18 Abril Kiruna 2024 
    '517': {0 : np.array([[0.951,  -0.612,	0.474	,0.459],
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

# Compute demodulation matrixes by inverting
demod_matrices = {}
for filt in mod_matrices:
    demod_matrices[filt] = {}
    for cam in mod_matrices[filt]:
        demod_matrices[filt][cam] = np.linalg.inv(mod_matrices[filt][cam])




# ------------------------------  CODE  ------------------------------------------ # 

"""
def demodulate(folder, dmod_matrices, filt, mode = 'single', wvl_index = 4):
    
    files = sorted(glob.glob(f"{folder}/*.img"))

    sorted_files = np.reshape(files, (3, 7, 4, 2))
    
    if filt == 'F517':
        filt_index = 0
    elif filt == 'F06':
        filt_index = 1
    elif filt == 'F02':
        filt_index = 2
    else:
        print(f"Filter must be: F517, F02 or F06")
 
    if mode == 'single':

        Npix = 1798 - 250
        data = np.zeros((2, 4, Npix, Npix))

        for ind, im in enumerate(sorted_files[filt_index, wvl_index, :, 0]):
            I, _ = read_Tumag(im)
            data[0, ind] = I[250:1798, 250:1798]
        for ind, im in enumerate(sorted_files[filt_index, wvl_index, :, 1]):
            I, _ = read_Tumag(im)
            data[1, ind] = I[250:1798, 250:1798]

        data[0] /= np.max(data[0])
        data[1] /= np.max(data[1])

        demod = np.zeros(np.shape(data))

        dm_cam1 = np.matmul(dmod_matrices[filt]["CAM1"], np.reshape(data[0], (4, Npix * Npix)))
        dm_cam2 = np.matmul(dmod_matrices[filt]["CAM2"], np.reshape(data[1], (4, Npix * Npix)))

        demod[0] = np.reshape(dm_cam1, (4, Npix, Npix))
        demod[1] = np.reshape(dm_cam2, (4, Npix, Npix))

        dual_beam = demod[0] * 0.5 + np.flip(demod[1], axis = 2) * 0.5

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

        dual_beam = demod[0] * 0.5 + np.flip(demod[1], axis = 2) * 0.5

    return data, demod, dual_beam
"""
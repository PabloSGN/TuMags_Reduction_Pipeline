# ---------------------------- DESCRIPTION --------------------------------------- #

# ------------------------------ IMPORTS ----------------------------------------- #

import os
import glob 
import time
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------ CONFIG ------------------------------------------ #

## Readout noise
fu_readout = [0.24, -0.24]
du_readout = [0.1, 0.1]
fv_readout = [0, 0]
dv_readout = [0.8, 0.8]

# ------------------------------  CODE  ------------------------------------------ # 

def create_frequency_mask(image_shape, fu, du, fv, dv):
    """
    Creates a binary mask to filter out specific frequencies, independent of image size.

    Parameters:
    - image_shape: Tuple (rows, cols) of the new image.
    - fu, fv: Center of the unwanted frequency (normalized).
    - du, dv: Bandwidth around the center frequency to remove.

    Returns:
    - Binary mask (same shape as image) with 0 at unwanted frequencies.
    """
    rows, cols = image_shape
    u = np.fft.fftfreq(rows)[:, None]  # Vertical frequency values (normalized)
    v = np.fft.fftfreq(cols)[None, :]  # Horizontal frequency values (normalized)

    # Create a mask initialized to 1 (keep all frequencies)
    mask = np.ones((rows, cols), dtype=np.float32)

    for i in range(len(fu)):
        # Apply the band-stop filter using a combined condition
        mask[(np.abs(u - fu[i]) < du[i]) & (np.abs(v - fv[i]) < dv[i])] = 0

    return mask

def filter_frecuencies(data, fu = fu_readout, du = du_readout, fv = fv_readout, dv = dv_readout, onelambda = False):
    
    """
    Filters an image by removing specified frequencies.

    Parameters:
    - data (np.array): data array dims (2, lambda, mods, Nx, Ny) 
    - fu (float, default : readout freq) : frequency to filter in the u dimension
    - fv (float, default : readout freq) : frequency to filter in the v dimension
    - du (float, default : readout freq) : margin to create the mask in the u dimension
    - dv (float, default : readout freq) : margin to create the mask in the v dimension
    - onelambda (Boolean, default : False) : Boolean in case only onelambda is passed

    Returns:
    - Filtered data.
    """
    if onelambda:
        data = data[:, np.newaxis] # To allow for only one lamdba.
    # Get shape for data
    shape = np.shape(data)
    nlambda = shape[1]
    nmods = shape[2]

    mask = create_frequency_mask(shape[-2:], fu, du, fv, dv)

    filtered = np.zeros(shape)

    print("Performing noise filtration ...")
    for lambd in range(nlambda):
        tic = time.time()
        print(f"Procesing wavelength: {lambd + 1} / {nlambda}")
        for mod in range(nmods):
            for cam in range(2):
                # Compute 2D Fourier Transform
                F = np.fft.fft2(data[cam, lambd, mod])
                F_shifted = np.fft.fftshift(F) # Center around 0
                F_filtered = F_shifted * mask # Mask slected frecuencies
                F_inverse_shifted = np.fft.ifftshift(F_filtered) # Return to original position
                filtered[cam, lambd, mod] = np.fft.ifft2(F_inverse_shifted).real  # Take real part
        tac = time.time()
        print(f"Time elapsed: {round(tac - tic, 3)} s.")

    print(f"Filtering process completed.\n")
   
    if onelambda:
        return filtered[:, 0]
    else:
        return filtered




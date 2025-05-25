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

def mirror_pad_5d(array, pad=100):
    """
    Pads a 5-dimensional array by mirroring the values along the last two dimensions (x and y).

    Parameters:
    -----------
    array : numpy.ndarray
        A 5-dimensional array to be padded. The dimensions are expected to follow the order (a, b, c, x, y).
    pad : int, optional
        The width of the padding to be applied to the last two dimensions (x and y). Default is 100.

    Returns:
    --------
    numpy.ndarray
        A 5-dimensional array with the last two dimensions padded using the 'reflect' mode.
    """
    # Pad only the last two dimensions (x and y)
    pad_width = [
        (0, 0),  # a
        (0, 0),  # b
        (0, 0),  # c
        (pad, pad),  # x
        (pad, pad),  # y
    ]
    return np.pad(array, pad_width=pad_width, mode='reflect')

def crop_5d(array, pad=100):
    """
    Crops a 5-dimensional array along the last two dimensions by removing a specified padding.

    Parameters:
        array (numpy.ndarray): A 5-dimensional array to be cropped. The expected shape is 
                               (batch, channel, depth, height, width).
        pad (int, optional): The number of pixels to remove from each side of the last two 
                             dimensions (height and width). Default is 100.

    Returns:
        numpy.ndarray: The cropped 5-dimensional array with the same number of dimensions, 
                       but reduced size in the last two dimensions.
    """
    return array[:, :, :, pad:-pad, pad:-pad]

def gaussian_band_mask(shape, axis='horizontal', center=None, width=10):
    """
    Generate a Gaussian band mask for an image.

    This function creates a 2D mask with a Gaussian band applied either 
    horizontally or vertically. The mask values range from 1 (no masking) 
    to values closer to 0 (masked regions) based on the Gaussian function.

    Parameters:
        shape (tuple): The shape of the mask as (rows, cols).
        axis (str, optional): The orientation of the Gaussian band. 
                              Options are 'horizontal' or 'vertical'. 
                              Default is 'horizontal'.
        center (int, optional): The center of the Gaussian band. If None, 
                                it defaults to the middle of the respective axis.
        width (float, optional): The standard deviation of the Gaussian 
                                 function, controlling the width of the band. 
                                 Default is 10.

    Returns:
        numpy.ndarray: A 2D array of the same shape as the input, containing 
                       the Gaussian band mask.
    """
    rows, cols = shape
    mask = np.ones(shape)

    if center is None:
        center = rows // 2 if axis == 'horizontal' else cols // 2

    if axis == 'horizontal':
        for i in range(rows):
            mask[i, :] *= 1 - np.exp(-((i - center) ** 2) / (2 * width**2))
    else:  # vertical band
        for j in range(cols):
            mask[:, j] *= 1 - np.exp(-((j - center) ** 2) / (2 * width**2))

    return mask

def find_peaks_in_spectrum(spectrum, verbose=False):
    """
    Identifies peaks in a 1D power spectrum while excluding a specified region 
    around the center of the spectrum.

    Parameters:
    -----------
    spectrum : numpy.ndarray
        A 1D array representing the power spectrum line cut.
    verbose : bool, optional
        If True, plots the spectrum, the modified spectrum, and the detected peaks 
        for visualization. Default is False.

    Returns:
    --------
    peaks : numpy.ndarray
        Indices of the detected peaks in the spectrum.

    Notes:
    ------
    - The function excludes a region of ±N points around the center of the spectrum 
        to avoid detecting peaks in that area.
    - Peaks are identified using the `scipy.signal.find_peaks` function with a 
        prominence threshold of 2 times the root mean square (RMS) of the modified spectrum 
        and a minimum distance of 300 points between peaks.
    - The excluded center region is replaced with the mean value of the spectrum 
        for peak detection purposes.
    - If `verbose` is True, a plot is displayed showing the original spectrum, 
        the modified spectrum, and the detected peaks.

    Example:
    --------
    >>> import numpy as np
    >>> spectrum = np.random.random(1000)  # Example spectrum
    >>> peaks = find_peaks_in_spectrum(spectrum, verbose=True)
    >>> print(peaks)
    """
    from scipy.signal import find_peaks

    # spectrum: 1D numpy array (your power spectrum line cut)
    # Exclude center: ignore ±N points around center index
    N = 250  # number of points to exclude around center
    center = len(spectrum) // 2
    # spectrum = gaussian_filter1d(spectrum, sigma=5)

    # Create a mask to ignore center
    mask = np.ones_like(spectrum, dtype=bool)
    mask[center - N:center + N] = False
    spectrum_ = np.where(mask, spectrum, np.mean(spectrum))  # Set center to mean

    # rms = np.sqrt(np.mean((spectrum_ - np.mean(spectrum_))**2))
    # threshold = np.mean(spectrum_) + 3 * rms
    rms = np.std(spectrum_)
    peaks, _ = find_peaks(spectrum_, prominence=2 * rms, distance=300)

    # peaks, _ = find_peaks(spectrum_, height=threshold)
    # peaks, _ = find_peaks(spectrum_, prominence=0.3)

    # Plot
    if verbose:
        plt.plot(spectrum)
        plt.plot(spectrum_)
        plt.plot(peaks, spectrum[peaks], "rx")
        plt.title("Detected Peaks (excluding center)")
        plt.show()
    return peaks

def filter_frecuencies(data, fu = fu_readout, du = du_readout, fv = fv_readout, dv = dv_readout,
    onelambda = False, verbose = False, band = None):
    """

    This function applies frequency filtering to a 5D data array. It supports two modes:
    1. Band-based filtering, where specific frequency bands are removed using Gaussian masks.
    2. General frequency filtering, where a mask is created based on specified frequency ranges.
    - data (np.array): Input data array with dimensions (2, lambda, mods, Nx, Ny).
    - fu (float, optional): Frequency to filter in the u dimension. Defaults to `fu_readout`.
    - fv (float, optional): Frequency to filter in the v dimension. Defaults to `fv_readout`.
    - du (float, optional): Margin to create the mask in the u dimension. Defaults to `du_readout`.
    - dv (float, optional): Margin to create the mask in the v dimension. Defaults to `dv_readout`.
    - onelambda (bool, optional): If True, processes only one lambda. Defaults to False.
    - verbose (bool, optional): If True, prints detailed information during processing. Defaults to False.
    - band (float, optional): Width of the Gaussian band for band-based filtering. If None, general frequency filtering is applied.
    - np.array: Filtered data array with the same shape as the input.
    Notes:
    - If `band` is provided, the function performs band-based filtering by identifying peaks in the power spectrum and applying Gaussian masks.
    - If `band` is not provided, the function creates a frequency mask based on the specified `fu`, `fv`, `du`, and `dv` values.
    - The function supports verbose mode for debugging and visualization of intermediate steps.
    - The input data is expected to have 5 dimensions: (2 cameras, lambda, mods, Nx, Ny).
    onelambda = False, verbose = False, band = None):
    
    Filters an image by removing specified frequencies.

    Parameters:
    - data (np.array): data array dims (2, lambda, mods, Nx, Ny) 
    - fu (float, default : readout freq) : frequency to filter in the u dimension
    - fv (float, default : readout freq) : frequency to filter in the v dimension
    - du (float, default : readout freq) : margin to create the mask in the u dimension
    - dv (float, default : readout freq) : margin to create the mask in the v dimension
    - onelambda (Boolean, default : False) : Boolean in case only onelambda is passed
    - verbose (Boolea, default : False) : Print info on terminal. 

    Returns:
    - Filtered data.
    """
    if onelambda:
        data = data[:, np.newaxis] # To allow for only one lamdba.
    # Get shape for data
    shape = np.shape(data)
    nlambda = shape[1]
    nmods = shape[2]

    if band:

        data = mirror_pad_5d(data, pad=100)
        shape = np.shape(data)

        fft = np.fft.fft2(data[0,0,0,:,:])
        fft_shifted = np.fft.fftshift(fft)
        power_spectrum_cam1 = np.abs(fft_shifted) ** 2

        fft = np.fft.fft2(data[1,0,0,:,:])
        fft_shifted = np.fft.fftshift(fft)
        power_spectrum_cam2 = np.abs(fft_shifted) ** 2
        peaks_cam1 = find_peaks_in_spectrum(np.mean(np.log10(power_spectrum_cam1),axis=1),verbose=verbose)
        peaks_cam2 = find_peaks_in_spectrum(np.mean(np.log10(power_spectrum_cam2),axis=1),verbose=verbose)
        mask_cam1 = np.ones((shape[-2], shape[-1]), dtype=np.float32)
        mask_cam2 = np.ones((shape[-2], shape[-1]), dtype=np.float32)
        for nm in peaks_cam1:
            mask_cam1 *= gaussian_band_mask(shape[-2:], axis='horizontal', center=nm, width=band)
        for nm in peaks_cam2:
            mask_cam2 *= gaussian_band_mask(shape[-2:], axis='horizontal', center=nm, width=band)
        if verbose:
            plt.imshow(mask_cam1, cmap='gray')
            plt.colorbar()
            plt.show()
            plt.imshow(mask_cam2, cmap='gray')
            plt.colorbar()
            plt.show()
        for lambd in range(nlambda):
            if verbose:
                print(f"Procesing filtration of wavelength: {lambd + 1} / {nlambda}")
            for mod in range(nmods):
                # Compute 2D Fourier Transform
                F = np.fft.fft2(data[0, lambd, mod])
                F_shifted = np.fft.fftshift(F) # Center around 0
                F_filtered = F_shifted * mask_cam1 # Mask slected frecuencies
                F_inverse_shifted = np.fft.ifftshift(F_filtered) # Return to original position
                data[0, lambd, mod] = np.fft.ifft2(F_inverse_shifted).real  # Take real part
                F = np.fft.fft2(data[1, lambd, mod])
                F_shifted = np.fft.fftshift(F) # Center around 0
                F_filtered = F_shifted * mask_cam2 # Mask slected frecuencies
                F_inverse_shifted = np.fft.ifftshift(F_filtered) # Return to original position
                data[1, lambd, mod] = np.fft.ifft2(F_inverse_shifted).real  # Take real part

        data = crop_5d(data, pad=100)

        tac = time.time()
        if verbose:
            print(f"Time elapsed: {round(tac - tic, 3)} s.")

        print(f"Filtering process completed.\n")
    
        if onelambda:
            return data[:, 0]
        else:
            return data

    else:

        mask = create_frequency_mask(shape[-2:], fu, du, fv, dv)

        filtered = np.zeros(shape)

        for lambd in range(nlambda):
            tic = time.time()
            if verbose:
                print(f"Procesing filtration of wavelength: {lambd + 1} / {nlambda}")
            for mod in range(nmods):
                for cam in range(2):
                    # Compute 2D Fourier Transform
                    F = np.fft.fft2(data[cam, lambd, mod])
                    F_shifted = np.fft.fftshift(F) # Center around 0
                    F_filtered = F_shifted * mask # Mask slected frecuencies
                    F_inverse_shifted = np.fft.ifftshift(F_filtered) # Return to original position
                    filtered[cam, lambd, mod] = np.fft.ifft2(F_inverse_shifted).real  # Take real part
            tac = time.time()
            if verbose:
                print(f"Time elapsed: {round(tac - tic, 3)} s.")

        print(f"Filtering process completed.\n")
    
        if onelambda:
            return filtered[:, 0]
        else:
            return filtered




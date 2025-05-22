# Alignment Module Documentation

## Description
This module contains functions related to the alignment of observation modes. It includes functionalities for subpixel realignment, field stop detection, image rotation, frequency filtering, and quadrant-based alignment.

Developed by: **Instituto de Astrofísica de Andalucía (IAA-CSIC)**

## Dependencies
The module requires the following libraries:

- `numpy`
- `time`
- `matplotlib`
- `scipy.fftpack`
- `scipy.ndimage`
- `pd_functions_v22` (custom module)
- `image_filtering` (custom module)

## Functions

### `dftreg(F, G, kappa)`
Calculates the shift between two images with subpixel accuracy using the method from Sicairos 2008.

#### Parameters:
- `F` (np.ndarray): FFT of the first image (without fftshift).
- `G` (np.ndarray): FFT of the second image (without fftshift).
- `kappa` (float): Inverse of subpixel precision (e.g., `kappa=20` results in `0.05` pixel precision).

#### Returns:
- `error` (float): Registration error.
- `row_shift` (float): Row shift value.
- `col_shift` (float): Column shift value.
- `Gshift` (np.ndarray): Shifted image in Fourier space.

---

### `dftups(M, n_out, kappa, roff, coff)`
Performs upsampled cross-correlation via matrix multiplication.

#### Parameters:
- `M` (np.ndarray): Input image in Fourier domain.
- `n_out` (int): Number of pixels in the upsampled DFT.
- `kappa` (float): Inverse of subpixel precision.
- `roff` (int): Row offset.
- `coff` (int): Column offset.

#### Returns:
- `np.ndarray`: Upsampled cross-correlation result.

---

### `FTpad(IM, Nout)`
Zero-pads an image in Fourier domain to upsample it.

#### Parameters:
- `IM` (np.ndarray): Input image in Fourier domain.
- `Nout` (int): Output size after padding.

#### Returns:
- `np.ndarray`: Padded image.

---

### `realign_subpixel(ima, accu=0.01, verbose=True, return_shift=False)`
Aligns a series of images with subpixel precision using the Sicairos method.

#### Parameters:
- `ima` (np.ndarray): 3D array of shape `(Nima, Nx, Ny)`.
- `accu` (float, default=0.01): Accuracy of alignment in pixel units.
- `verbose` (bool, default=True): Print alignment process.
- `return_shift` (bool, default=False): Whether to return shift values.

#### Returns:
- `np.ndarray`: Aligned image array.
- `(Optional) row_shifts, col_shifts` (lists): Shifts applied to each image.

---

### `find_fieldstop(cam1=None, verbose=False, plot_flag=False, margin=10)`
Finds the field stop in an image.

#### Parameters:
- `cam1` (np.ndarray): Single image.
- `verbose` (bool, default=False): Print processing information.
- `plot_flag` (bool, default=False): Plot field stop detection.
- `margin` (int, default=10): Margin pixels from detected field-stop.

#### Returns:
- `np.ndarray`: Field stop coordinates.

---

### `rotate_camera2(data, theta=0.065, onelambda=False)`
Rotates camera 2 in an observation mode dataset.

#### Parameters:
- `data` (np.ndarray): Observation mode array `(Ncams, Nlambda, Nmods, Nx, Ny)`.
- `theta` (float, default=0.065): Rotation angle.
- `onelambda` (bool, default=False): Set to `True` if only one lambda is used.

#### Returns:
- `np.ndarray`: Rotated dataset.

---

### `filter_and_rotate(data, theta=0.0655, verbose=False, filterflag=True, zkes=np.zeros(21))`
Applies filtering and rotation to camera 2.

#### Parameters:
- `data` (np.ndarray): Observation mode array.
- `theta` (float, default=0.0655): Rotation angle.
- `verbose` (bool, default=False): Print process.
- `filterflag` (bool, default=True): Apply Fourier filtration.
- `zkes` (np.ndarray, default=`np.zeros(21)`): Zernike coefficients for filtration.

#### Returns:
- `np.ndarray`: Filtered and rotated dataset.

---

### `align_obsmode(data, acc=0.01, verbose=False, theta=0.0655, filterflag=True, onelambda=False, zkes=np.zeros(21), returnshifts=False)`
Filters, rotates, and aligns an observation mode.

#### Parameters:
- `data` (np.ndarray): Observation mode array.
- `acc` (float, default=0.01): Alignment accuracy.
- `theta` (float, default=0.0655): Rotation angle.
- `verbose` (bool, default=False): Print process.
- `filterflag` (bool, default=True): Apply Fourier filtration.
- `onelambda` (bool, default=False): Handle single lambda datasets.
- `zkes` (np.ndarray, default=`np.zeros(21)`): Zernike coefficients.
- `returnshifts` (bool, default=False): Return shift values.

#### Returns:
- `np.ndarray`: Aligned dataset.
- `(Optional) shifts` (list): Shift values.

---

### `reshape_into_16_quadrants(images, nlambda, nmods)`
Reshapes an observation mode into 16 quadrants (4x4 grid).

#### Parameters:
- `images` (np.ndarray): Observation mode array.
- `nlambda` (int): Number of wavelengths.
- `nmods` (int): Number of modulations.

#### Returns:
- `np.ndarray`: Reshaped dataset.

---

### `align_quadrants(data, acc=0.01, verbose=False)`
Aligns quadrants separately in an observation mode.

#### Parameters:
- `data` (np.ndarray): Observation mode array.
- `acc` (float, default=0.01): Alignment accuracy.
- `verbose` (bool, default=False): Print process.

#### Returns:
- `np.ndarray`: Aligned dataset.
- `shifts` (list): Shift values applied per quadrant.

---

## Usage
To use the alignment functions, import the module and call the desired function. Example:

```python
from alignment import align_obs_mode

aligned_images = align_obs_mode(obs_mode_data, accu=0.01)
```

---


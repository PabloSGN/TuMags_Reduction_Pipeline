# Demodulation Module Documentation

## Description
This module contains functions for performing demodulation on observation mode data. It includes predefined demodulation matrices and methods for processing full images as well as individual quadrants.

## Dependencies
The module requires the following libraries:

- `numpy`

## Predefined Matrices

The module includes predefined modulation and demodulation matrices for different filters (`517`, `525.02`, `525.06`). These matrices are used to compute the inverse demodulation matrices:

- `mod_matrices`: Standard modulation matrices.
- `mod_matrices_david`: Alternative modulation matrices.
- `demod_matrices`: Inverted matrices for demodulation.
- `demod_matrices_david`: Alternative inverted matrices.

## Functions

### `demodulate(data, filt, dmod_matrices=demod_matrices_david, onelambda=False, BothCams=False, verbose=False)`

Performs demodulation of observation mode data.

#### Parameters:
- `data` (np.ndarray): Observation mode array `(Ncams, Nlambda, Nmods, Nx, Ny)`.
- `filt` (str): Filter for demodulation (`517`, `525.02`, or `525.06`).
- `dmod_matrices` (dict, default=`demod_matrices_david`): Demodulation matrices to use.
- `onelambda` (bool, default=False): Whether the data contains only one wavelength.
- `BothCams` (bool, default=False): Whether to return individual camera outputs in addition to the combined output.
- `verbose` (bool, default=False): Whether to print debug information.

#### Returns:
- `np.ndarray`: Dual-beam demodulated data `(Nlambda, Nmods, Nx, Ny)`.
- `(Optional) np.ndarray`: Individual camera demodulated data `(Ncams, Nlambda, Nmods, Nx, Ny)` if `BothCams=True`.

---

### `demodulate_quadrants(data, nlambda, nmods, filt, nquads=16, Np_quad=354)`

Performs demodulation on observation mode data that has been split into quadrants.

#### Parameters:
- `data` (np.ndarray): Observation mode array `(Ncams, Nlambda, Nmods, Nquads, Nx, Ny)`.
- `nlambda` (int): Number of wavelengths.
- `nmods` (int): Number of modulations.
- `filt` (str): Filter for demodulation (`517`, `525.02`, or `525.06`).
- `nquads` (int, default=16): Number of quadrants.
- `Np_quad` (int, default=354): Number of pixels per quadrant.

#### Returns:
- `np.ndarray`: Dual-beam demodulated data `(Nlambda, Nmods, Nquads, Nx, Ny)`.
- `np.ndarray`: Individual camera demodulated data `(Ncams, Nlambda, Nmods, Nquads, Nx, Ny)`.

## Usage

To use the demodulation functions, import the module and call the desired function. Example:

```python
from demodulation import demodulate

dual_beam_output = demodulate(data, filt='517')
```

---
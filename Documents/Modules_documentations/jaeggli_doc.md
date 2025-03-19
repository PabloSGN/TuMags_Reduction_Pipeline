Here is the documentation for the `jaeggli.py` file, following the format and style of the provided `image_filtering_doc.md` template:

---

# Jaeggli Module Documentation

## Description
This module contains functions to correct for crosstalk using the method described in Jaeggli et al. 2022. The module implements diattenuation and retarder modeling to minimize the correlation between Stokes parameters (I, Q, U, V) in weakly polarized regions of spectral lines. The method is based on the paper [Jaeggli et al. 2022](https://doi.org/10.3847/1538-4357/ac6506) and the associated [GitHub repository](https://github.com/sajaeggli/adhoc_xtalk).

## Dependencies
The module requires the following libraries:

- `numpy`
- `matplotlib`
- `glob`
- `scipy.optimize.minimize`
- `pandas`

## Functions

### `polmodel1(D, theta, chi)`

Computes the diattenuation Mueller matrix based on the given parameters.

#### Parameters:
- `D` (float): Diattenuation value.
- `theta` (float): Polar angle.
- `chi` (float): Azimuthal angle.

#### Returns:
- `np.ndarray`: Diattenuation Mueller matrix.

---

### `fitfunc1(param, stokesin)`

Computes the merit function for the diattenuation model based on the input Stokes vector and parameters.

#### Parameters:
- `param` (tuple): Tuple containing `(D, theta, chi)`.
- `stokesin` (np.ndarray): Input Stokes vector.

#### Returns:
- `float`: Merit function value.

---

### `minimize_for_model1(iMM, bs)`

Computes the merit function for a given inverse Mueller matrix and Stokes vector.

#### Parameters:
- `iMM` (np.ndarray): Inverse Mueller matrix.
- `bs` (np.ndarray): Input Stokes vector.

#### Returns:
- `float`: Merit function value.

---

### `polmodel2(theta, delta)`

Computes the retarder Mueller matrix based on the given parameters.

#### Parameters:
- `theta` (float): Polar angle.
- `delta` (float): Retardance.

#### Returns:
- `np.ndarray`: Retarder Mueller matrix.

---

### `fitfunc2(fitangles, stokesin)`

Computes the merit function for the retarder model based on the input Stokes vector and parameters.

#### Parameters:
- `fitangles` (tuple): Tuple containing `(theta, delta)`.
- `stokesin` (np.ndarray): Input Stokes vector.

#### Returns:
- `float`: Merit function value.

---

### `minimize_for_model2(iMM, bs)`

Computes the merit function for a given inverse Mueller matrix and Stokes vector.

#### Parameters:
- `iMM` (np.ndarray): Inverse Mueller matrix.
- `bs` (np.ndarray): Input Stokes vector.

#### Returns:
- `float`: Merit function value.

---

### `fit_mueller_matrix(data, pthresh=0.02, norm=False, region=[200, 1200, 200, 1200], last_wvl=None, plots=False)`

Fits the best diattenuation Mueller matrix to minimize the correlation between Stokes I and Q, U, V in weakly polarized regions.

#### Parameters:
- `data` (np.ndarray): 4D array with Stokes parameters (wavelength, Stokes, x, y).
- `pthresh` (float, default=0.02): Stokes V polarization threshold to determine weak/strong polarization regions.
- `norm` (bool, default=False): Whether to normalize Stokes components to the median of Stokes I at the first wavelength.
- `region` (list, default=[200, 1200, 200, 1200]): Region of the image to consider for the fit (x0, xf, y0, yf).
- `last_wvl` (int or None, default=None): Last wavelength to consider for the merit function.
- `plots` (bool, default=False): Whether to plot the fractional polarization map.

#### Returns:
- `datarest` (np.ndarray): 4D array with corrected Stokes parameters.
- `MM1a` (np.ndarray): Diattenuation Mueller matrix.

---

## Usage

To use the module, import the desired functions and call them with the appropriate parameters. Example:

```python
from jaeggli import fit_mueller_matrix

# Fit the Mueller matrix and correct the data
corrected_data, mueller_matrix = fit_mueller_matrix(data, pthresh=0.02, plots=True)
```
---

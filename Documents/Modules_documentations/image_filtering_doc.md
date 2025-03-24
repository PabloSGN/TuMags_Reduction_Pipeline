# Image Filtering Module Documentation

## Description
This module provides functions to filter out specific frequencies from images using Fourier transforms. It includes predefined frequency settings for readout noise filtering and functions for generating frequency masks and filtering images.

## Dependencies
The module requires the following libraries:

- `os`
- `glob`
- `time`
- `numpy`
- `matplotlib`

## Predefined Frequency Settings

The module defines readout noise parameters for filtering:
- `fu_readout`: Center frequencies in the u-dimension.
- `du_readout`: Bandwidths in the u-dimension.
- `fv_readout`: Center frequencies in the v-dimension.
- `dv_readout`: Bandwidths in the v-dimension.

## Functions

### `create_frequency_mask(image_shape, fu, du, fv, dv)`

Creates a binary mask to filter out specific frequencies, independent of image size.

#### Parameters:
- `image_shape` (tuple): Shape `(rows, cols)` of the image.
- `fu` (list): Center frequencies in the u-dimension (normalized).
- `du` (list): Bandwidths around `fu` to remove.
- `fv` (list): Center frequencies in the v-dimension (normalized).
- `dv` (list): Bandwidths around `fv` to remove.

#### Returns:
- `np.ndarray`: Binary mask of the same shape as the image, with `0` at unwanted frequencies.

---

### `filter_frecuencies(data, fu=fu_readout, du=du_readout, fv=fv_readout, dv=dv_readout, onelambda=False)`

Applies a frequency filter to remove specified frequencies from the input data.

#### Parameters:
- `data` (np.ndarray): Input array `(2, Nlambda, Nmods, Nx, Ny)`.
- `fu` (list, default=`fu_readout`): Frequencies to filter in the u-dimension.
- `fv` (list, default=`fv_readout`): Frequencies to filter in the v-dimension.
- `du` (list, default=`du_readout`): Margins for filtering in the u-dimension.
- `dv` (list, default=`dv_readout`): Margins for filtering in the v-dimension.
- `onelambda` (bool, default=False): Whether the data contains only one wavelength.

#### Returns:
- `np.ndarray`: Filtered data with unwanted frequencies removed.

## Usage

To use the filtering functions, import the module and call the desired function. Example:

```python
from image_filtering import filter_frecuencies

filtered_data = filter_frecuencies(data)
```

---



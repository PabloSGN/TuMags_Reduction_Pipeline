Here is the documentation for the `master_flatfield.py` file, following the format and style of the provided `image_filtering_doc.md` template:

---

# Master Flatfield Module Documentation

## Description
This module provides functionality to compute the master flat field from a set of flat-field observations for a single observation mode. It includes options for normalization and prefilter removal.

## Dependencies
The module requires the following libraries:

- `numpy`
- `time`
- `config` (custom module for configuration settings)
- `image_handler` (custom module for image handling)
- `prefilter_removal` (custom module for prefilter removal)

## Functions

### `compute_master_flat_field(flat_fields_paths, dc, lambda_repeat=4, verbose=False, norm_method="avg", remove_prefilter=False, pref_model=None, volts=None)`

Computes the master flat field from a set of flat-field observations. The function normalizes the flat fields and optionally removes the prefilter contribution.

#### Parameters:
- `flat_fields_paths` (list): List of paths to the flat-field images.
- `dc` (np.ndarray): Dark current array.
- `lambda_repeat` (int, default=4): Lambda repeat parameter of the observation.
- `verbose` (bool, default=False): Whether to print progress and details.
- `norm_method` (str, default="avg"): Normalization method. Options are `"avg"` (average of all modulations) or `"mod"` (each modulation separately).
- `remove_prefilter` (bool, default=False): Whether to remove the prefilter contribution from the flat fields.
- `pref_model` (optional): Prefilter model file required if `remove_prefilter=True`.
- `volts` (str or None, default=None): Set to `"read"` to use read voltages for prefilter removal. If `None`, fixed voltages are used.

#### Returns:
- `ff_data` (np.ndarray): Array containing the master flat field (cams, Nlambda, Nmods, Nx, Ny).
- `ff_info` (dict): Dictionary containing all information about the flat-field observation.

---

## Usage

To use the module, import the `compute_master_flat_field` function and call it with the appropriate parameters. Example:

```python
from master_flatfield import compute_master_flat_field

# Compute the master flat field
flat_field_data, flat_field_info = compute_master_flat_field(flat_fields_paths, dc, verbose=True, remove_prefilter=True, pref_model="prefilter_model.pkl")
```

---

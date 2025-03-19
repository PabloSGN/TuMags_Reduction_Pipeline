Here is the documentation for the `master_dark.py` file, following the format and style of the provided `image_filtering_doc.md` template:

---

# Master Dark Module Documentation

## Description
This module provides functionality to compute the dark current from a set of dark observations. It averages all observations and returns the dark current per accumulation for both cameras.

## Dependencies
The module requires the following libraries:

- `numpy`
- `utils` (custom module for reading Tumag files)
- `time`

## Functions

### `compute_master_darks(dark_paths, verbose=False)`

Computes the dark current from a set of dark observations. The function averages all dark images and returns the dark current per accumulation for both cameras.

#### Parameters:
- `dark_paths` (list): List of paths to the dark current images.
- `verbose` (bool, default=False): Whether to print progress and details.

#### Returns:
- `np.ndarray`: Array of the dark current for both cameras, normalized per accumulation.

---

## Usage

To use the module, import the `compute_master_darks` function and call it with the appropriate parameters. Example:

```python
from master_dark import compute_master_darks

# Compute the master dark current
dark_current = compute_master_darks(dark_paths, verbose=True)
```
---


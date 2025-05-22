# Image Handler Module Documentation

## Description
This module provides functionality to handle and process `.img` files. It includes classes for reading image headers and metadata, managing observations, and processing images for various observation modes. The module is designed to work with data from the Tumag instrument, including dark current correction, image binning, and organization of images by observation counters and modes.

## Dependencies
The module requires the following libraries:

- `os`
- `time`
- `matplotlib.pyplot`
- `numpy`
- `pandas`
- `datetime`
- `utils` (custom module for reading Tumag files)
- `config` (custom module for configuration settings)

## Classes and Functions

### `raw_header` Class

This class creates a header object for images from raw values. It includes transfer functions to convert counts to volts for various components (HVPS, LCVR) and stores metadata about the image.

#### Methods:

- **`hvps_commanded_2_volts(counts: int, sign: int) -> float`**  
  Converts HVPS commanded counts to volts.  
  **Parameters:**  
  - `counts` (int): The count value to convert.  
  - `sign` (int): The sign of the voltage.  
  **Returns:**  
  - `float`: The converted voltage value.

- **`hvps_read_2_volts(counts: int) -> float`**  
  Converts HVPS read counts to volts.  
  **Parameters:**  
  - `counts` (int): The count value to convert.  
  **Returns:**  
  - `float`: The converted voltage value.

- **`lcvr1_counts_2_volts(counts: int) -> float`**  
  Converts LCVR1 counts to volts.  
  **Parameters:**  
  - `counts` (int): The count value to convert.  
  **Returns:**  
  - `float`: The converted voltage value.

- **`lcvr2_counts_2_volts(counts: int) -> float`**  
  Converts LCVR2 counts to volts.  
  **Parameters:**  
  - `counts` (int): The count value to convert.  
  **Returns:**  
  - `float`: The converted voltage value.

- **`__init__(camID: int, om: int, nAcc: int, Roix: int, Roiy: int, Roix_offset: int, Roiy_offset: int, oc: int, fw1: int, fw2: int, hvps_counts: int, hvps_sign: int, lcvr1_counts: int, lcvr2_counts: int, hvps_read_counts: int, lcvr1_read_counts: int, lcvr2_read_counts: int, image_name: str)`**  
  Initializes the header object with metadata from the image.  
  **Parameters:**  
  - `camID` (int): Camera ID.  
  - `om` (int): Observation mode index.  
  - `nAcc` (int): Number of accumulations.  
  - `Roix`, `Roiy` (int): ROI dimensions.  
  - `Roix_offset`, `Roiy_offset` (int): ROI offsets.  
  - `oc` (int): Observation counter.  
  - `fw1`, `fw2` (int): Filter wheel positions.  
  - `hvps_counts`, `hvps_sign` (int): HVPS commanded counts and sign.  
  - `lcvr1_counts`, `lcvr2_counts` (int): LCVR counts.  
  - `hvps_read_counts`, `lcvr1_read_counts`, `lcvr2_read_counts` (int): Read counts for HVPS and LCVRs.  
  - `image_name` (str): Name of the image file.

- **`get_info()`**  
  Returns the metadata dictionary.  
  **Returns:**  
  - `dict`: Dictionary containing image metadata.

---

### `read(image_path: str)`

Reads a single `.img` file and returns the image data and header information.

#### Parameters:
- `image_path` (str): Path to the `.img` file.

#### Returns:
- `np.ndarray`: Image data.
- `dict`: Header information.

---

### `nominal_observation` Class

This class processes observation modes, organizing images by wavelength and modulation, and correcting for dark current.

#### Methods:

- **`__init__(om, images_path, dc)`**  
  Initializes the observation mode object.  
  **Parameters:**  
  - `om` (int): Observation mode.  
  - `images_path` (list): List of paths to the images.  
  - `dc` (np.ndarray): Dark current array.

- **`get_info()`**  
  Returns the observation mode information.  
  **Returns:**  
  - `dict`: Observation mode metadata.

- **`get_data()`**  
  Returns the processed image data.  
  **Returns:**  
  - `np.ndarray`: Processed image data.

---

### `nominal_flat` Class

This class processes flat-field observation modes, averaging multiple repetitions and correcting for dark current.

#### Methods:

- **`__init__(om, images_path, nreps, dc, lambda_repeat=4, verbose=False)`**  
  Initializes the flat-field mode object.  
  **Parameters:**  
  - `om` (int): Observation mode.  
  - `images_path` (list): List of paths to the images.  
  - `nreps` (int): Number of repetitions.  
  - `dc` (np.ndarray): Dark current array.  
  - `lambda_repeat` (int, default=4): Number of wavelength repetitions.  
  - `verbose` (bool, default=False): Whether to print progress.

- **`get_info()`**  
  Returns the flat-field mode information.  
  **Returns:**  
  - `dict`: Flat-field mode metadata.

- **`get_data()`**  
  Returns the processed flat-field data.  
  **Returns:**  
  - `np.ndarray`: Processed flat-field data.

---

### `get_images_paths(queries)`

Provides the paths of the images associated with the queries. Queries should be in the format `"DXX-start-end"`, where `DXX` is the observation day.

#### Parameters:
- `queries` (str or list): Query or list of queries.

#### Returns:
- `list`: List of image paths.

---

### `read_ID(image_index, plotflag=False, verbose=False, header=False, binning=False)`

Reads an image with a given ID and optionally plots, bins, or prints header information.

#### Parameters:
- `image_index` (str): Image ID in the format `"DXX-YYYY"`.
- `plotflag` (bool, default=False): Whether to plot the image.
- `verbose` (bool, default=False): Whether to print verbose information.
- `header` (bool, default=False): Whether to print the header.
- `binning` (bool, default=False): Whether to bin the image.

#### Returns:
- `np.ndarray`: Image data.
- `dict`: Header information.

---

### `separate_ocs(paths, verbose=True, flat_fieldmode=False)`

Separates images into different observation counters (OCs).

#### Parameters:
- `paths` (list): List of image paths.
- `verbose` (bool, default=True): Whether to print progress.
- `flat_fieldmode` (bool, default=False): Whether to enable flat-field mode.

#### Returns:
- `dict`: Dictionary containing organized images by OCs.

---

### `get_time_from_filename(filename)`

Extracts the timestamp from the filename and returns it as a `datetime` object.

#### Parameters:
- `filename` (str): Name of the image file.

#### Returns:
- `datetime`: Timestamp extracted from the filename.

---

### `obs_mode_separator(paths, verbose=False)`

Organizes images of an observation mode based on header properties.

#### Parameters:
- `paths` (list): List of image paths.
- `verbose` (bool, default=False): Whether to print progress.

#### Returns:
- `dict`: Dictionary containing organized images by observation mode.

---

### `check_timestamps(paths, verbose=True)`

Checks the intervals between images and generates a plot showing the intervals.

#### Parameters:
- `paths` (list): List of image paths.
- `verbose` (bool, default=True): Whether to print progress.

---

### `snapshot_processing(paths, dc, verbose=True)`

Processes snapshot images (HC timelines).

#### Parameters:
- `paths` (list): List of image paths.
- `dc` (np.ndarray): Dark current array.
- `verbose` (bool, default=True): Whether to print progress.

#### Returns:
- `np.ndarray`: Processed snapshot data.

---

### `polarizers_parser(paths, filt, dc)`

Parses micropolarizer observations.

#### Parameters:
- `paths` (list): List of image paths.
- `filt` (str): Filter to parse (`525.02`, `517`, `525.06`, or `all`).
- `dc` (np.ndarray): Dark current array.

#### Returns:
- `np.ndarray`: Processed micropolarizer data.

---

## Usage

To use the module, import the desired functions or classes and call them with the appropriate parameters. Example:

```python
from image_handler import read, nominal_observation

# Read an image
image_data, header_info = read("path/to/image.img")

# Process an observation mode
observation = nominal_observation(om="1", images_path=get_images_paths("D10-1000-2000"), dc=dark_current_array)
```
---

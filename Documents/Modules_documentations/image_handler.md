
---

## **Description**
This module provides functionality to handle and process `.img` files, specifically for reduced observations in astrophysics. It includes classes for reading image headers and metadata, performing calculations, and managing observational data.

**Author:**  
Pablo Santamarina Guerrero  
Instituto de Astrofísica de Andalucía (IAA-CSIC)  
Email: pablosantamarinag@gmail.com  

---

## **Imports**

```python
import time
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from utils import read_Tumag
import config as cf
```

- **time**: Used for tracking time in processing tasks.
- **os**: Handles directory and file path operations.
- **matplotlib.pyplot**: Used for displaying images.
- **numpy**: Supports numerical computations, especially array manipulations.
- **pandas**: Helps with data manipulation and reading `.csv` files.
- **read_Tumag**: A function from the custom `utils` module for reading `.img` files.
- **config**: Custom configuration module `config` as `cf`, which includes observation modes and filter wheel configurations.

---

## **Classes and Functions**

### **`class raw_header`**
Handles the conversion of raw header values to physical quantities and creates a header dictionary for each image.

#### **Methods**
1. **`__init__(...)`**  
   Initializes the `raw_header` object with the image metadata and converts counts to volts for various components.
   
   **Parameters:**
   - `camID`: Camera ID.
   - `om`: Observation Mode index.
   - `nAcc`: Accumulation count.
   - `Roix`, `Roiy`: Size of image.
   - `Roix_offset`, `Roiy_offset`: ROI offsets.
   - `oc`: Observation Counter.
   - `fw1`, `fw2`: Filter wheel positions.
   - `hvps_counts`, `hvps_sign`: High Voltage Power Supply counts and sign.
   - `lcvr1_counts`, `lcvr2_counts`: LCVR counts.
   - `hvps_read_counts`: HVPS read values in counts.
   - `lcvr1_read_counts`, `lcvr2_read_counts`: LCVR real counts.

2. **`hvps_commanded_2_volts(counts, sign)`**  
   Converts HVPS commanded counts to volts.

3. **`hvps_read_2_volts(counts)`**  
   Converts HVPS read counts to volts.

4. **`lcvr1_counts_2_volts(counts)`**  
   Converts LCVR1 counts to volts.

5. **`lcvr2_counts_2_volts(counts)`**  
   Converts LCVR2 counts to volts.

6. **`get_info()`**  
   Returns the header information dictionary with converted voltages and original metadata.

---

### **`def read(image_path: str)`**
Reads a single `.img` file, extracts its header, and creates a `raw_header` object with relevant image metadata.  

**Parameters:**
- `image_path`: Path to the `.img` file.

**Returns:**
- `img`: The image data.
- `head.get_info()`: Header information extracted from the image.

---

### **`class nominal_observation`**
Processes multiple `.img` files for a given observation mode and organizes the data into arrays.

#### **Methods**
1. **`__init__(om, images_path, dc)`**  
   Initializes the `nominal_observation` class by reading images and headers, subtracting dark current, and organizing data into arrays for each camera and wavelength.

   **Parameters:**
   - `om`: Observation mode.
   - `images_path`: List of image file paths.
   - `dc`: Dark current to subtract from the images.

2. **`get_info()`**  
   Returns the header information and observation metadata.

3. **`get_data()`**  
   Returns the data array of the images.

---

### **`class nominal_flat`**
Similar to `nominal_observation`, but processes flat field observations, typically over multiple repetitions for calibration purposes.

#### **Methods**
1. **`__init__(om, images_path, nreps, dc, lambda_repeat=4, verbose=False)`**  
   Processes images across repetitions and wavelength repeats, subtracts dark current, and stores results in arrays.

   **Parameters:**
   - `om`: Observation mode.
   - `images_path`: List of image paths.
   - `nreps`: Number of repetitions.
   - `dc`: Dark current to subtract.
   - `lambda_repeat`: Number of wavelength repeats (default is 4).
   - `verbose`: If `True`, prints additional progress information.

2. **`get_info()`**  
   Returns the header information and observation metadata.

3. **`get_data()`**  
   Returns the data array of the images.

---

### **`def get_images_paths(queries)`**
Generates a list of image paths based on queries that specify the day and range of indices.

**Parameters:**
- `queries`: A list or single query in the format `"DXX-start-end"`, where `DXX` is the day (e.g., "D09") and `start`, `end` are the image index range.

**Returns:**
- `selection`: List of file paths matching the query.

---

### **`def read_ID(image_index, plotflag=False, verbose=False, header=False)`**
Reads an image by its index and optionally plots the image, prints header information, or both.

**Parameters:**
- `image_index`: The identifier of the image, including the day and index.
- `plotflag`: If `True`, displays the image.
- `verbose`: If `True`, prints observation counter and mode.
- `header`: If `True`, prints the header information.

**Returns:**
- `I`: The image data.
- `H`: The image header.

---

### **`def separate_ocs(paths, verbose=True)`**
Groups images by Observation Counter (OC) and verifies the completeness of each observation set.

**Parameters:**
- `paths`: List of image file paths.
- `verbose`: If `True`, prints progress and results.

**Returns:**
- `OCs`: Dictionary where each OC has information about the observation mode, expected number of images, and completeness status.

---

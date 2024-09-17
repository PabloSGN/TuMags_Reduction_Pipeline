Here is the documentation for the provided Python module:

# **Module Documentation**

## **Description**
This module provides a function to compute the master flat field from a set of flat-field observations for a single observation mode. The function reads a series of flat-field images, corrects them for dark current, normalizes them, and returns the corrected flat-field data.

## **Imports**
- **`numpy (np)`**: Used for numerical computations, specifically handling arrays and matrices.
- **`time`**: Used for timing the execution of the function.
- **`utils.read_Tumag`**: A custom function (not used directly in this function) for reading `.img` files.
- **`config (cf)`**: Contains observation mode configuration data such as the number of wavelengths and modulations for each observation mode.
- **`image_handler (ih)`**: A custom module that handles images, including processing flat-field modes..

## **Configuration**
The configuration settings are implicitly handled within the function through the `config` module, specifically for retrieving observation mode configurations.

## **Function: `compute_master_flat_field`**

### **Parameters:**
- **`flat_fields_paths`**: A list of file paths for flat-field images. These paths should correspond to images taken in the same observation mode.
- **`dc`**: The dark current data that will be used to correct the flat-field images. This is scaled by the number of accumulations.
- **`lambda_repeat`**: An integer (default = 4) that defines how many times each wavelength is repeated during the observation.
- **`verbose`**: A boolean flag (`True/False`). If `True`, the function prints detailed output, including the number of flat images, observation mode, repetitions, wavelengths, and modulations.

### **Returns:**
- **`norm_ff`**: A numpy array containing the normalized flat-field data.
- **`flat_obs.get_info()`**: Metadata information related to the flat-field observations.

### **Description:**
1. **Initialization and Verbose Output**: The function starts by recording the time and printing initial information (if `verbose` is enabled). It displays the number of flat-field images to be processed.
  
2. **Header Information**: The function reads the header of the first flat-field image to retrieve important metadata, such as:
   - Observation mode (`ObservationMode`).
   - Number of wavelengths (`N_wls`) and modulations (`N_mods`) for the observation mode, fetched from the `config` module.

3. **Repetitions Check**: It checks whether the number of provided flat-field images is divisible by the expected number of images based on the number of wavelengths, modulations, and `lambda_repeat`. If not, an error is raised indicating incomplete observations.

4. **Flat-field Image Reading**: The function reads the flat-field images using the `ih.nominal_flat` function, applying a dark current correction (scaled by the number of accumulations).

5. **Normalization of Flat-fields**: After retrieving the flat-field data, the function normalizes it by dividing each image by its mean value (excluding the outermost 300 pixels) for both cameras.

6. **Completion**: The function prints the total time taken for the computation (if `verbose` is enabled) and returns the normalized flat-field data (`norm_ff`) and metadata (`flat_obs.get_info()`).

### **Example Usage:**
```python
flat_fields_paths = image_handler.get_images_paths("D10-100-700")
dark_current = compute_master_darks(dark_paths)
norm_ff, flat_info = compute_master_flat_field(flat_fields_paths, dark_current, verbose=True)
```

### **Output Example:**
```
Computing flats.
------------------
N flats: 40
Observation Mode: OM1
Nº of repetitions: 10
Nº of wavelengths: 5
Nº of Modulations: 2
Flat-fields computed in 3.456 s.
```
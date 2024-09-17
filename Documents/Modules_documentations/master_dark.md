Here is the documentation for the provided Python module:

# **Module Documentation**

## **Description**
This module provides a function to compute the dark current from a set of observations. The function reads a series of images, averages them, and returns the dark current per accumulation to allow scaling for different observations. 

The function also supports optional verbosity to print detailed information about the process.

## **Imports**
- **`numpy (np)`**: Used for numerical computations, specifically for handling arrays and matrices.
- **`utils.read_Tumag`**: Custom function to read `.img` files, assumed to return an image and header metadata.
- **`time`**: Used for timing the execution of the function.


## **Configuration**
There is no explicit configuration section in this module. Configuration is handled implicitly by the input parameters.

## **Function: `compute_master_darks`**

### **Parameters:**
- **`dark_paths`**: A list of file paths for dark current images. These paths should correspond to two cameras (`"_0_"` for camera 1, and `"_1_"` for camera 2).
- **`verbose`**: A boolean flag (`True/False`). If `True`, the function will print detailed output, including the number of dark images and the number of accumulations.

### **Returns:**
- **`dark_current`**: A 3D numpy array of shape `(2, x, y)`, where:
  - The first index (`0` or `1`) represents the camera (1 or 2).
  - The `x` and `y` dimensions represent the image size.
  - The values represent the averaged dark current per pixel, scaled by the number of accumulations.

### **Description:**
1. **Image Size Determination**: The function reads the first image from `dark_paths` to determine the size of the image (`x` by `y` dimensions).
2. **Dark Current Initialization**: An empty array is initialized to hold the dark current data for both cameras (`dark_current[0]` for camera 1 and `dark_current[1]` for camera 2).
3. **Dark Image Separation**: The list of dark image paths is split into two groups—one for camera 1 (`darks_cam_1`) and one for camera 2 (`darks_cam_2`), based on the presence of the `_0_` or `_1_` identifier in the filename.
4. **Dark Current Calculation**: 
    - For each image in `darks_cam_1`, the image data is read and added to the corresponding index of `dark_current[0]`.
    - For each image in `darks_cam_2`, the image is read, flipped horizontally, and added to `dark_current[1]`.
5. **Averaging**: The total dark current for each camera is averaged by dividing by the number of dark images for that camera.
6. **Per Accumulation Scaling**: The resulting dark current is divided by the number of accumulations (retrieved from the header of the first image).
7. **Execution Time**: The function prints the total time taken for computation if `verbose` is enabled.

### **Example Usage:**
```python
dark_paths = image_handler.get_images_paths("D10-100-199")
dark_current = compute_master_darks(dark_paths, verbose=True)
```
OR:

```python
dark_paths = glob.glob("/path/to/dark/images/*.img")
dark_current = compute_master_darks(dark_paths, verbose=True)
```

### **Output Example:**
```
Computing darks.
------------------
N darks for cam1 : 10
N darks for cam2 : 10
N accumulations : 100
Dark current computed in 2.345 s.
```
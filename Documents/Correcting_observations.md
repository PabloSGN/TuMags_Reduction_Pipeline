# Correcting observations

The standard procedure of dtaa correction is the following:
1. Compute the **dark current** (DC).
2. Compute the **flat-fields** (FFs).
3. Compute the field-stop from the FFs and align the two cameras. 
4. Process the observation mode and align the cameras.
5. Correct the FFs. 
6. Demodulate. 
7. Correct X-talk. 

Let's go over these steps...

## Identifiying images. 

First identify the images corresponding to each calibration block using the [Tumag's Logsheet](https://docs.google.com/spreadsheets/d/1RJ5KIgxMN6B-1xDe9gRfoTh1L_uTDMbZw0ajLK6S0so/edit?usp=sharing)

Let's say we want to process the minimum success observations (Comm 1 + Comm 2).

The images are identified through the following indexes: 
```python
flats_obs1_index = "D10-2740-4339"
flats_obs2_index = "D10-4340-5619"
darks_index = "D10-5620-5719"
pinholes_index = "D10-6188-6193"
obs_index = "D10-304-2439"
```

All routines of the pipeline require the paths to the images. There are two ways of working with this: 

1. <ins>Working in the IAA's server.</ins>

If working on the IAA's server, the paths can be obtained directly with the indexes and the routine [get_images_paths](../image_handler.py#L230) from the module [image_handler.py](../image_handler.py) : 

```python
dark_paths = ih.get_images_paths(darks_index)
ff_obs1_paths = ih.get_images_paths(flats_obs1_index)
ff_obs2_paths = ih.get_images_paths(flats_obs2_index)
pinholes_paths = ih.get_images_paths(pinholes_index)
obs_images = ih.get_images_paths(obs_index)
```
2. <ins>Working locally.</ins>

If working on a local machine you will need to download the data, and generate the paths of your machine. <ins>**Remember to provide them sorted!** </ins>

If stored in separate folders this can be achieved with the [glob](https://docs.python.org/3/library/glob.html) module by: 
```python
darks_paths = sorted(glob.glob("path/to/darks/directory/*")
```
There are other ways but remember to sort them out. 

## Compute the **dark current**

To generate the dark current use the function [compute_master_dark](../master_dark.py#L18) from the module [master_dark.py](../master_dark.py):

```python
dc = compute_master_darks(dark_paths, verbose = True)  
```
The compute_master_darks only requires the images paths as an input. 

The dark current (variable dc) is an array of shape (2 x 2016 x 2016) - (camera, size x,size y). 

The dark current is given per acumulation. 

<ins>Note:</ins> The verbose entry for all functions is used to print info on the terminal. Set it to False or don't use it if no info is wanted.

## Compute the **flat-fields**

To generate the flat-field use the function [compute_master_flat_field](../master_flatfield.py#L15) from the module [master_flatfield.py](../master_flatfield.py):

```python
ff_data, ff_info = compute_master_flat_field(ff_paths, dc = dc, verbose = True)
```

The function compute_master_flat_fields requires the paths to the flats observation and the dark current. Only flat-fields of one observing mode can be provided for each execution of the function. 

<ins>Note:</ins> Currently, if one observing mode is incomplete (there is some image missing) the routine will fail. This will be upgraded...

The function returns two variables:
1. ff_data : A numpy array of shape (2 x Nlambda x Nmods x 2016 x 2016) with the normalized flat fields.
2. ff_info : A dictionary containing information of the observing mode and images headers.

## Compute alignment between cameras. 

In order to combine the two images from the two cameras, the images must be aligned. 



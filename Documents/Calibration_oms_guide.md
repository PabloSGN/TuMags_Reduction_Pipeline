# Calibration observing modes guide

Disclaimer: As for now there are no specific functions to process linear polarizers, micropolarizers and prefilters. 

The TuMag's calibration modes are:
1. [Dark current](#dark-current)
2. [Flat-fields](#flat-fields)
3. [Phase diversity](#phase-diversity)
4. [Linear polarizers](#linear-polarizers)
5. [Micropolarizers](#Micropolarizers)
6. [Straylight](#straylight)
7. [Prefilter scans](#prefilter-scans)
8. [Pinholes](#pinholes)

## Dark current

Dark current observations are obscured observations composed of 100 images (50 per camera). 

Total number of images: 100
Observation Mode: Spectral Scan

<ins>Processing images example:<ins>

Using the function: "compute_master_dark" from module [master_Dark](../master_dark.py#L18):
```python
dc = master_dark.compute_master_dark("D10-5620-5719", verbose = True) 
```
This commands processes all images from the observations and returns the average of all images per camera and per acumulation and stores them in a numpy array of shape: (N cameras (2) x sizex x sizey)

## Flat fields

Flat field modes are the same as nominal observation modes with multiple repetitions per wavelength (Nlambda = 4) and whole scans (Nreps). 

The number of images depends on the observation mode and the number of repetitions. 

<ins>Processing images example:<ins>
To generate the flat-field use the function [compute_master_flat_field](../master_flatfield.py#L15) from the module [master_flatfield.py](../master_flatfield.py):

```python
ff_data, ff_info = compute_master_flat_field(ff_paths, dc = dc, verbose = True)
```

The function compute_master_flat_fields requires the paths to the flats observation and the dark current. Only flat-fields of one observing mode can be provided for each execution of the function. 

<ins>Note:</ins> Currently, if one observing mode is incomplete (there is some image missing) the routine will fail. This will be upgraded...

The function returns two variables:
1. ff_data : A numpy array of shape (2 x Nlambda x Nmods x 2016 x 2016) with the normalized flat fields.
2. ff_info : A dictionary containing information of the observing mode and images headers.

## Phase diversity

Phase diveristy calibrtionsets consist on 4 different observation modes: 2 for each observed filter (depens on timeline) and one foces and one defocused. 

Each observation mode consists on 40 images per camera (start of the mission) or 32 images (changed in the middle of campaign, from D12-47445 onwards) per mode. 

Total number of images: 320 (start) or  256 (from D12-47445)
Observation Mode: PD calibration

<ins>Processing images example:<ins>

Using the function: "process_pd_observation" from module [phase_diversity](../phase_diversity.py#L75):
```python
pd_data_fe = pd.process_pd_observation("D10-6194-6513", filt = 0, verbose = True) 
```
This commands reads the data from the first filter (filt = 0) and stores it in a numpy array of shape: (Ncameras (2) x Focus vs defocus (2) x N images x size x x size y).

## Linear polarizers

Observations employed using the linear polarizer. The three filters are usually measured for the vectorial mode and only two for thelongitudinal one. 

1. **Vectorial modes:**

4 images (4 modulations) per camera and per filter (3 prefilters). Total images: 24.
Observation Mode: Polarimetric calibration

2. **Longitudinal mode:**
2 images (2 modulations) per camera and per filter (2 prefilters). Total images: 8
Observation Mode: Polarimetric calibration

## Micropolarizers

Observations using the micropolarizers target. Always in vectorial configuration. 

4 images (4 modulations) per camera and per filter (3 prefilters). Total images: 24.
Observation Mode: Polarimetric calibration

## Straylight

Observation mode 1 employing the straylight target. Nimages : 80. 
Observation Mode: 1

<ins>Processing images example:<ins>

Since it is a nominal observation mode 1, the images can be processed with the *nominal_observation* function from the [image handler](../image_handler.py) module:
```python
straylight = ih.nominal_observation("1", "D10-7264-7343", dark_current) 
```
## Prefilter scans

Scans in voltage to measure the prefilters shape. Can be manual (launched manually to retireve thumbnails during the operations) or programmed in timelines. 
Observation Mode: Spectral calibration

1. **Programmed**

Number of images: 35 images per camera. 2 repetitions per filter (2 filters). Total: 280 images.
Observation Mode: Polarimetric calibration

## Pinholes

Pinhole observations are observations using the pinhole target of the filter wheel. 

Images: 2 per filter ( Nimages: 6)
Observation Mode: Spectral calibration
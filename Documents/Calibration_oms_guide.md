# Calibration observing modes guide

Disclaimer: As for now there are no specific functions to process linear polarizers, micropolarizers and prefilters. 

The TuMag's calibration modes are:
1. [Phase diversity](#phase-diversity)
2. [Linear polarizers](#linear-polarizers)
3. [Micropolarizers](#Micropolarizers)
4. [Straylight](#straylight)
5. [Prefilter scans](#prefilter-scans)

## Phase Diversity

Phase diveristy calibrtionsets consist on 4 different observation modes: 2 for each observed filter (depens on timeline) and one foces and one defocused. 

Each observation mode consists on 40 images per camera (start of the mission) or 32 images (changed in the middle of campaign, from D12-47445 onwards) per mode. 

Total number of images: 320 (start) or  256 (from D12-47445)

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

2. **Longitudinal mode:**
2 images (2 modulations) per camera and per filter (2 prefilters). Total images: 8


## Micropolarizers

Observations using the micropolarizers target. Always in vectorial configuration. 

4 images (4 modulations) per camera and per filter (3 prefilters). Total images: 24.

## Straylight

Observation mode 1 employing the straylight target. Nimages : 80. 

<ins>Processing images example:<ins>

Since it is a nominal observation mode 1, the images can be processed with the *nominal_observation* function from the [image handler](../image_handler.py) module:
```python
straylight = ih.nominal_observation("1", "D10-7264-7343", dark_current) 
```
## Prefilter scans

Scans in voltage to measure the prefilters shape. Can be manual (launched manually to retireve thumbnails during the operations) or programmed in timelines. 

1. **Programmed**

Number of images: 35 images per camera. 2 repetitions per filter (2 filters). Total: 280 images.
# TuMag's Pipeline Version tracker.

## v1.0
First stable version.

### Main features

#### Flat-field computation 
    - Flat-fields are normalized to the average of all modulations of each wavelength.  
    - Flat-fields are rescaled to account for the pre-filter when correcting the data based on a fit to the prefilter scan calibration mode.  
#### Alignment
    - Camera's 2 rotated by fixed value. 
    - Alignment routine permorms a filtering of the readout frequencies only (**not noise filtering**, this is only performed when using wavefront reconstruction)
    - Alignment based on the sicarios method.
    - For each wavelength, all modulations and cameras are aligned with respect to the first modulation of camera 0.
#### Cross - Talk correction
    - Standard cross-talk correction performed with jaeggli module (FJB approach). 

## v1.1
Changed output of demodulation and alignment.

### Demodulation
Demodualtion no longer returns both cameras unless specified. 

### Alignment
Alignment no longer returns shifts unless specified. 

# Planned improvements:

 - Implement prefilter fitting into flat-field calculation function seamlessly. 
 - Blueshift and calculation and correction. 
 - Overall improvement of flat-field computation.
 - Inclusion of automatic rotation angle for camera 2. 


## v1.2

### Config
- Added Pipeline version for automatic labeling. 

### Demodulation
- Added camera intenisty compensation for the dual beam. 

### fits_files_handling
- Added function to generate fits and headers. 

### main_processing
- Added function to use for the dtaa processing and fits generation. 

### phase_diversity
 - Added function to apply wavefront reconstruction. 

### reduction_example_generate_fits.py 
- Example of use of the main processing module.

### alignment.py
- Automatic rotation angle selection from interpolation.
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

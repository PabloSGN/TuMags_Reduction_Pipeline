# TuMag's reduction from scratch 

This is a guide of TuMag's reduction pipeline from scrath, directly from the .img files. For a guide of the reduction using the fits files see [Woriking with fits files](./Working_with_fits_files.md)


## Required data. 

TuMag's reduction pipeline is thought to process an individual [Observing Mode](../README.md#observation-modes-om) at a time. Fort each observing mode, the following data has to be identified:
 1. Closest dark current observations. 
 2. Closest and corresponding flat-field.
 3. Images corresponding to the observing mode to reduce. 

Some other data can also be employed for optional steps.
Optional:
 4. Closest and corresponding pre-filter scan (for prefilter removal)
 5. Associated Zernike coefficients (for wavelength reconstruction)

To learn how to identify the .img files see [Working with TuMag's data](./Working_with_Tumags_data.md)


# Example of calibration image IDs for Timeline EMF_1 AR2 block (Sunrise_ID : 02_EMEF_TM_01)
```python
dc_paths = image_handler.get_images_paths("D10-6844-6943")
ff_mode_1_paths = image_handler.get_images_paths("D10-7398-9637")
```

## Reduction steps.

TuMag's reduction pipeline consist on the following steps: 
 1.  [Compute the dark-current](#compute-the-dark-current). 
 2. [Compute the master flat-field.](#compute-the-master-flat-field) 
 3. [Process the observing mode](#process-the-observing-mode)
 4. [Aply the flat-field correction.](#aply-the-flat-field-correction)
 5. [Filter and align.](#filter-and-align) 
 6. [Demodulation](#demodulation)
 7. [Cross-talk correction.](#cross-talk-correction) 

Let's go over these steps in detail. Or jump to the [Summary](#summary)

### Compute the dark-current.

The dark-current is computed with the [compute_master_darks](../master_dark.py#L16) funtion from the [master_dark](../master_dark.py) module.

**Main arguments**
 1. dark_paths : list containing the paths to the dark current observations. 

**Optional argument**
 2. verbose (default False): boolean to print info on terminal. 

**Output**
 1. dc : dark-current observation (np.array -> dims : (2, Nx, Ny)) 

**EXAMPLE**
```python
dc = compute_master_darks(dark_paths = dc_paths, verbose = True)
```

### Compute the master flat-field.

The master flat-field is computed with the [compute_master_flat_field](../master_flatfield.py#L20) funtion from the [master_flatfield](../master_flatfield.py) module.

**Main arguments**
 1. flat_fields_paths : list containing the paths to the flat field observations. 
 2. dc : dark_current as computed in [previous step](#compute-the-dark-current). 

**Optional argument**
 3. lambda_repeat (default 4): Number of repetitions per wavelength. 
 4. (default False): boolean to print info on terminal.  
 5. norm_method (default "avg") : Noormalization method, can be "avg" to normalize the the average of all moduulations or "mod" to normalize each modulation by its mean. 
 6. remove_prefilter (default  False) : activate prefilter removal from flat-fields. 
 7. pref_model (default None) : Pickle file containing the prefilter model. 
 8. volts (default None) : if "read" read voltages are used in the prefilter removal. 

**Output**
 1. norm_ff : master flat-field (np.array -> dims : (2, Nlambda, Nmods, Nx, Ny)) 
 2. header_ff : Dictionary containing info of the flat-field. 

**EXAMPLE**
```python
ff_mode1, ff_mode1_info = master_flatfield.compute_master_flat_field(ff_mode_1_paths, dc = dc, verbose = True)
```

For a more detailed explanaition of the prefilter removal parameters see [Prefilter removal guide](./Fit_Prefilter_guide.md)

### Process the observing mode.

The observing modes are handled with the [nominal_observation](../image_handler.py#L109) class from the [image_handler](../image_handler.py) module.

**Main arguments**
 1. Observing Mode : String identifiyng the [observing mode](../README.md#observation-modes-om) to reduce.
 2. images_paths : List containing the paths to the observing mode images.
 3. dc : dark_current as computed in [step one](#compute-the-dark-current).

**Output**
 1. om : Observing mode class containg both data and header (np.array -> dims : (2, Nx, Ny)) 

The data can be extracted with the .get_data() method, and the header with the .get_info() method. 

Images belonging to the observing mode have to be previously identified (check [Working with TuMag's data](./Working_with_Tumags_data.md) for more info). 

**EXAMPLE**
```python
om = ih.nominal_observation(observing_mode_id, observing_mode_images_paths, dc = dc) 
om_data = om.get_data()
om_info = om.get_info()
```
### Aply the flat-field correction.

The flat-field correction can be computed with the [correct_observation](../master_flatfield.py#L115) function from the [master_flatfield](../master_flatfield.py) module.

**Main arguments**
 1. Observing Mode data : Numpy array contaiining the observing mode data the (np.array -> dims : (2, Nlambda, Nmods, Nx, Ny) / (2,Nmods, Nx, Ny) if one_lambda = True)) 
 2. Flat-field data : Numpy array contaiining the coputed flat-field (np.array -> dims : (2, Nlambda, Nmods, Nx, Ny)

**Optional arguments**
 3. onelambda (default : False) : Activate the possibility of using only one wavelength. 

**Output**
 1. corr_om : Corrected observing mode (np.array -> dims : (2, Nlambda, Nmods, Nx, Ny) / (2,Nmods, Nx, Ny) if one_lambda = True)

**EXAMPLE**
```python
om_corr = master_flatfield.correct_observation(data = om_data, ff = ff_mode1)
```

### Filter and align

Before demodulating the data, the observing mode has to be filtered to remove some spurious signals and aligned to be able to combine both cameras and different modulations. 

Both tasks are performed with the [align_obsmode](../alignment.py#L291) function from the [alignment](../alignment.py) module.

**Main arguments**
 1. Observing Mode data : Numpy array containing the observing mode data the (np.array -> dims : (2, Nlambda, Nmods, Nx, Ny) / (2,Nmods, Nx, Ny) if one_lambda = True)) 

**Optional arguments**
 2. acc (Default = 0.01) : Accuracy of the alignment method.
 3. verbose (default : False): boolean to print info on terminal. 
 4. theta (default : 0.0655) : Angle of camera's 2 rotation. 
 5. filterflag (default : True) : Activate frequency filter. 
 6. onelambda (default : False) : Activate the possibility of using only one wavelength. 
 7. returnshifts (default : False) : Activate rewturn of shifts of images.  

**Output**
 1. aligned : Aligned and filtered observing mode (np.array -> dims : (2, Nlambda, Nmods, Nx, Ny) / (2,Nmods, Nx, Ny) if one_lambda = True)
 2. shifts (if returnshifts = True) : shifts applied to each frame. 

**EXAMPLE**
```python
aligned = alignment.align_obs_mode(om_corr, verbose = True)
```

### Demodulation

After alignment, the data can be demodulated with the [demodulate](../demodulation.py#96) function from the [demodulation](../demodulation.py) module.

**Main arguments**
 1. data : Numpy array containing the observing mode data (np.array -> dims : (2, Nlambda, Nmods, Nx, Ny) / (2,Nmods, Nx, Ny) if one_lambda = True)) 
 2. filt = pre-filter corresponding to the observing mode (517 / 525.02 / 525.06, can be obtained from obs mode info.) 

**Optional arguments**
 2. dmod_matrices : demodulation matrixes.
 3. onelambda (default : False) : Activate the possibility of using only one wavelength. 
 4. BothCams (default : False) : Activate rewturn of both cameras demodulation.
 5. verbose (default : False): boolean to print info on terminal.  

**Output**
 1. dual_beam : Dual-beamed demodulated stokes parameters (np.array -> dims : (Nlambda, Nmods, Nx, Ny) / (Nmods, Nx, Ny) if one_lambda = True)
 2. demod (if BothCams = True) : demodulated stokes for each camera. (np.array -> dims : (2, Nlambda, Nmods, Nx, Ny) / (2,Nmods, Nx, Ny) if one_lambda = True)
    

**EXAMPLE**
```python
stokes = demodulation.demodulate(data = aligned, filt = om_info["line"] )
```
### Cross-talk correction

After demodulation, the cross-talk correction can be performed with the [fit_mueller_matrix](../xtalk_jaeggli.py#L131) function from the [xtalk_jaeggli](../xtalk_jaeggli.py) module.


**Main arguments**
 1. data : Numpy array containing the dual-beamed demodulated stokes parameters (np.array -> dims : (Nlambda, Nmods, Nx, Ny)) 

**Output**
 1. corrected_stokes : Corrected stokes parameters (np.array -> dims : (Nlambda, Nmods, Nx, Ny))
 2. MM1a : Fitted muller matrix of the method

**EXAMPLE**
```python
xtalk_corr = xtalk_jaeggli.fit_mueller_matrix(stokes)
```

## Summary 

The complete code for reducing a single observing mode using the approach of the Observation counters pickle file described in [Working with TuMag's data](./Working_with_Tumags_data.md).

```python

# Observing IDs
dc_paths = image_handler.get_images_paths("D10-6844-6943")
ff_mode_1_paths = image_handler.get_images_paths("D10-7398-9637")
OCs_pickle = "path/to/ocs/picklefile.pickle
oc_number = 176 # For example


with open(OCs_pickle, "rb") as file:
    OCS = pickle.load(file) 

# Compute darks
dc = compute_master_darks(dark_paths = dc_paths, verbose = True)

# Compute flats
ff_mode1, ff_mode1_info = master_flatfield.compute_master_flat_field(ff_mode_1_paths, dc = dc, verbose = True)

# Process observing mode
om = ih.nominal_observation(OCS[oc_number]["OM"], OCS[oc_number]["ims"], dc = dc) 
om_data = om.get_data()
om_info = om.get_info()

# Correct flat-fields
om_corr = master_flatfield.correct_observation(data = om_data, ff = ff_mode1)

# Align obs mode
aligned = alignment.align_obs_mode(om_corr, verbose = True)

# Demodulate
stokes = demodulation.demodulate(data = aligned, filt = om_info["line"] )

#Correct xtalk
xtalk_corr = xtalk_jaeggli.fit_mueller_matrix(stokes)
```
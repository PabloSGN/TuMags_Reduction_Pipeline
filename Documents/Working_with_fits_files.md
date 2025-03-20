# Working with Fits files. 

When working with TuMag's fits files, the first step is to identify the reduction level of the provided file. 

## Reduction levels

Fits files at different stages of the reduction process will be provided, identified with the **reduction level**. Additionally, data will be provided with and without wavefront reconstruction. The details of each level and the steps that have been applied to the data are summarized in this table: 

| Reduction Level | Flats + Darks | Wavefront reconstruction | Alignment and filtering | Stokes Parameters |
|:-:|:-:|:-:|:-:|:-:|
| **0.5** | :white_check_mark: | :x: | :x: | :x: |
| **0.6** | :white_check_mark: | :white_check_mark: | :x: | :x: |
| **0.8** | :white_check_mark: | :x: | :white_check_mark: | :x: |
| **0.9** | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x: |
| **1.0** | :white_check_mark: | :x: | :white_check_mark: | :white_check_mark: |
| **1.1** | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: |

## Working with TuMag's Fits files. 

TuMag's fits file names have the following structure: 

**{ Sunrise_ID } _ TM _ { BlockIndex } _ { Line }{ Observing Mode } _ { Nº Wavelengths } _ { Tstart } _ LV _ { reduction level} _ v{pipeline version}.fits**

 - Sunrise_ID : Common ID that identifies the timeline for the three instruments.
 - BlockIndex : Index (00, 01, 02, etc) that identifies the different observing blocks within the timeline for TuMag. 
 - Line : Mg or Fe
 - Observing Mode : [Observing Modes](../README.md#observation-modes-om) of TuMag
 - Nº Wavelengths : Number of wavelengths
 - Tstart : Time of the first image of the observing mode. Format : ddmmYYYYTHH:MM:SS 
 - reduction level : [Reduction level](./TuMags_reduction_guide.md#reduction-levels)
 - Pipeline version : Version of the pipeline generated to use the fits. 

Example: 01_QSUN_TM_00_Fe2.02_8_10072024T133713_LV_1.0_v1.0.fits

### Fits headers

All fits files have in the header the relevant information of each observing modes namely:
 - Spectral line
 - Observing mode.
 - Wavelengths of the scan.
 - Volts of the scan
 - Number of accumulations
 - Reduction level 
 - Pipeline version 

For reduction levels higher than 0.5 some fields are added. 

#### Zernikes 
If wavefront reconstruction has been applied to the data, the zernike coefficients are included in the main header. 

#### Shifts
For any fits of level 0.8 and higher, the shifts for all wavelengths, modulations and cameras are included in the first extension of the header under the name "SHIFTS"

#### Muller Matrix
For any file of level 1.0 and higher the fitted muller matrix of the jaeggli module is also included in the second header extension under the name "Fitted_Muller_Matrix"

## Reduction from fits files. 

Data is stored in the fits file in an array of shape : (Ncams, Nlambda, Nmods, Nx, Ny) for any reduction level lower than 1.0, and (Nlambda, Nstokes, Nx, Ny) for higher levels. To read them, use astropy's fits module:

```python
from astropy.io import fits
data = fits.getdata(fits_file)
header = fits.getheader(fits_file)
```

Depending on the initial reduction level of the file, some steps are required before the stokes-parameters are ready to be used. Since all data has at been dark-current corrected and flat-fielded, these steps are: 
 1. [Filter and align.](#filter-and-align---for-files-with-reduction-level-lower-than-08) 
 2. [Demodulation](#demodulation)
 3. [Cross-talk correction.](#cross-talk-correction) 


### Filter and align - For files with reduction level lower than 0.8

Before demodulating the data, the observing mode has to be filtered to remove some spurious signals and aligned to be able to combine both cameras and different modulations. 

Both tasks are performed with the [align_obsmode](../alignment.py#L291) function from the [alignment](../alignment.py) module.

**Main arguments**
 1. Observing Mode data : Numpy array containing the observing mode data.

**Optional arguments**
 - acc (Default = 0.01) : Accuracy of the alignment method.
 - verbose (default : False): boolean to print info on terminal. 
 - theta (default : 0.0655) : Angle of camera's 2 rotation. 
 - filterflag (default : True) : Activate frequency filter. 
 - onelambda (default : False) : Activate the possibility of using only one wavelength. 
 - returnshifts (default : False) : Activate rewturn of shifts of images.  

**Output**
 1. aligned : Aligned and filtered observing mode (np.array -> dims : (2, Nlambda, Nmods, Nx, Ny) / (2,Nmods, Nx, Ny) if one_lambda = True)
 2. shifts (if returnshifts = True) : shifts applied to each frame. 

**EXAMPLE**
```python
aligned = alignment.align_obs_mode(data, verbose = True)
```

### Demodulation - For files with reduction level lower than 1.0

After alignment, the data can be demodulated with the [demodulate](../demodulation.py#96) function from the [demodulation](../demodulation.py) module.

**Main arguments**
 1. data : Numpy array containing the aligned observing mode data. 
 2. filt = pre-filter corresponding to the observing mode (517 / 525.02 / 525.06, can be obtained from header.) 

**Optional arguments**
 - dmod_matrices : demodulation matrixes.
 - onelambda (default : False) : Activate the possibility of using only one wavelength. 
 - BothCams (default : False) : Activate rewturn of both cameras demodulation.
 - verbose (default : False): boolean to print info on terminal.  

**Output**
 1. dual_beam : Dual-beamed demodulated stokes parameters (np.array -> dims : (Nlambda, Nmods, Nx, Ny) / (Nmods, Nx, Ny) if one_lambda = True)
 2. demod (if BothCams = True) : demodulated stokes for each camera. (np.array -> dims : (2, Nlambda, Nmods, Nx, Ny) / (2,Nmods, Nx, Ny) if one_lambda = True)
    

**EXAMPLE**
```python
stokes = demodulation.demodulate(data = aligned, filt = header["FW2"])
```
### Cross-talk correction - For files with reduction level lower than 1.0 

After demodulation, the cross-talk correction can be performed with the [fit_mueller_matrix](../xtalk_jaeggli.py#L131) function from the [xtalk_jaeggli](../xtalk_jaeggli.py) module.


**Main arguments**
 1. data : Numpy array containing the dual-beamed demodulated stokes parameters. 

**Output**
 1. corrected_stokes : Corrected stokes parameters (np.array -> dims : (Nlambda, Nmods, Nx, Ny))
 2. MM1a : Fitted muller matrix of the method

**EXAMPLE**
```python
xtalk_corr = xtalk_jaeggli.fit_mueller_matrix(stokes)
```

## Summary 

The complete code for reducing a single observing mode is:

```python
# Read the fits file
data = fits.getdata(fits_file)
header = fits.getheader(fits_file)

# Align obs mode
aligned = alignment.align_obs_mode(data, verbose = True) # For reduction level < 0.8

# Demodulate
stokes = demodulation.demodulate(data = aligned, filt = header["FW2"] )  # For reduction level < 1.0

#Correct xtalk
xtalk_corr = xtalk_jaeggli.fit_mueller_matrix(stokes) # For reduction level < 1.0
```

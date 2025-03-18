# TuMag's reduction step by step guide. 


## Reduction levels

Fits files at different stages of the reduction process will be provided, identified with the **reduction level**. Additionally, data will be provided with and without wavefront reconstruction. The details of each level and the steps that have been applied to the data are summarized in this table: 

| Reduction Level | Flats + Darks | Wavefront reconstruction | Alignment and filtering | Stokes Parameters |
|----------|----------|----------|----------|----------|
| **0.5** | :white_check_mark: | :x: | :x: | :x: |
| **0.6** | :white_check_mark: | :white_check_mark: | :x: | :x: |
| **0.8** | :white_check_mark: | :x: | :white_check_mark: | :x: |
| **0.9** | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x: |
| **1.0** | :white_check_mark: | :x: | :white_check_mark: | :white_check_mark: |
| **1.1** | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: |

## Working with TuMag's Fits files. 

TuMag's fits file names have the following structure: 

**{ SunriseID }_TM_{ BlockIndex } _ { Line }{ Observing Mode }_{ Nº Wavelengths } _ { Tstart } _ LV_{ reduction level} _ v{pipeline version}.fits**

 - SunriseID : Common ID that identifies the timeline for the three instruments.
 - BlockIndex : Index (00, 01, 02, etc) that identifies the different observing blocks within the timeline for TuMag. 
 - Line : Mg or Fe
 - Observing Mode : [Observing Modes](../README.md#observation-modes-om) of TuMag
 - Nº Wavelengths : Number of wavelengths
 - Tstart : Time of the first image of the observing mode. Format : ddmmYYYYTHH:MM:SS 
 - reduction level : [Reduction level](./TuMags_reduction_guide.md#reduction-levels)
 - Pipeline version : Version of the pipeline generated to use the fits. 
 
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
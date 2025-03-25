# TuMag's pipeline

This is the github for TuMag's reduction pipeline. 

## Overview of the main modules for the reduction:
 1. [**alignment.py**](./alignment.py) : Module with functions to find the fieldstop and align the frames of observing modes. [Doc](Documents/Modules_documentations/alignment_doc.md)
 2. [**demodulation.py**](./demodulation.py) : Module with the functions to demodulate the data. [Doc](Documents/Modules_documentations/demodulation_doc.md)
 3. [**image_filtering.py**](./image_filtering.py) : Module to filter spurious frequencies from the image. [Doc](Documents/Modules_documentations/image_filtering_doc.md)
 4. [**image_handler.py**](./image_handler.py) : Module to manage .img files. [Doc](Documents/Modules_documentations/image_handler_doc.md)
 5. [**master_dark.py**](./master_dark.py) : Module to process dark current observations. [Doc](Documents/Modules_documentations/master_dark_doc.md)
 6. [**master_flatfield.py**](./master_flatfield.py) : Module to process flat-field observations. [Doc](Documents/Modules_documentations/master_flat_field_doc.md)
 7. [**xtalk_jaeggli.py**](./xtalk_jaeggli.py) : Module to compute the cross-talk from observations. [Doc](Documents/Modules_documentations/xtalk_jaeggli_doc.md)

## Some useful guides
- [Installation](Documents/Installation.md)
- [Working with TuMag's data](Documents/Working_with_Tumags_data.md)
- [Correcting observations from fits files](Documents/Working_with_fits_files.md)
- [Correcting observations from fits files - Jupyter](./TuMags_reduction_example.ipynb)
- [Correcting observations from scratch](Documents/TuMags_Reduction_from_scratch.md)
- [Prefilter Removal](Documents/Fit_Prefilter_guide.md)

## TuMag's data paroperties 

When working with TuMag's data some key concepts should be known. 

### Timelines

All observations (scientific and calibrations) are recorded and grouped in their corresponding Timeline in the [Tumag's Logsheet](https://docs.google.com/spreadsheets/d/1RJ5KIgxMN6B-1xDe9gRfoTh1L_uTDMbZw0ajLK6S0so/edit?usp=sharing).  There are copies of the logsheet provided as a [pdf](Documents/TuMags%20Logsheet%20-%20Timelines%20detailed.pdf) and [csv](Documents/TuMags%20Logsheet%20-%20Timelines%20detailed.csv).

TuMag observations were stuctured through  **timelines**, a list of pre-configured commands where calibration and scientific observations were already programmed. 

For a guide describing TuMag's commands for each timeline see [Julián's webpage](https://www.uv.es/jublanro/TuMag_timeline_reference.html). 

Go to [Tumag's official data website](https://www.uv.es/jublanro/tumag_data_test.html) for already reduced observations. 

A summary of the observations, the corresponding timelines and pointing information can be found in [Observation summary](Documents/Sunrise_summary_observations.pdf).

<ins>Important:</ins> When working with the data have in mind that the timelines were paused and stopped, some commands were run manually, others aborted, etc ... 

### Observation Modes (OM)

Tumag observations are operated through pre-configured **observation modes**. 

<ins>Science obsering modes:</ins> 

| **Observing Mode** | **Filter** | **Nº wavelengths** | **Modulation scheme** | 
|:--------:|:--------:|:--------:|:--------:|
| 0p    | 517    | 12  | Vectorial   | 
| 0s    | 517    | 12  | No modulation   |
| 1    | 517    | 10 | Vectorial  | 
| 2.02    | 525.02   | 8  | Vectorial  |
| 2.06    | 525.06  | 8  | Vectorial  |
| 3.02    | 525.02   | 5  | Longitudinal  |
| 3.06    | 525.06  | 5  | Longitudinal  |
| 4    | 517  | 3  | Vectorial  |
| 5.02    | 525.02   | 3  | Vectorial  |
| 5.06    | 525.06  | 3  | Vectorial  |

<ins>Calibration obsering modes:</ins> 

Calibration observing modes are more flexible than science observing modes and can be configured to match the observations. For a detailed description of each calibration OM see [Calibration OMs guide](Documents/Calibration_oms_guide.md)

### Observation Counters

The **observation counter** (OC) is a field in the [images' header](Documents/Image_header.md) that identifies the observation mode corresponding to the image. All images of a specific observation mode share the same OC. The counter goes from 0 to 255 and then cycles.

Usually observations are identified with this field in [TuMag's Logbook](Documents/TuMags%20Logsheet%20-%20Timelines%20detailed.csv).

## Contributing

IAA-SPG team and Julian.

## Versions

For a log of the version and the changes see: [versions Log](./version_log.md)

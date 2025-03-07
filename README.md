# TuMag's pipeline

Tools for reading and processing raw data from TuMag.

Disclaimer: Very early still...

## Useful guides
- [Installation](Documents/Installation.md)
- [Working with TuMag's data](Documents/Working_with_Tumags_data.md)
- [Correcting observations](Documents/Correcting_observations.md)
- [Module guide](Documents/Module_guide.md)

## Overview of TuMag's observations

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

Calibration observing modes are more flexible than science observing modes nd can be confugured to match the observations. For a detailed description of each calibration OM see [Calibration OMs guide](Documents/Calibration_oms_guide.md)

### Observation Counters

The **observation counter** (OC) is a field in the [images' header](Documents/Image_header.md) that identifies the observation mode corresponding to the image. All images of a specific observation mode share the same OC. The counter goes from 0 to 255 and then cycles.

Usually observations are identified with this field in [TuMag's Logbook](Documents/TuMagCompass.csv).

## Contributing

SPG team.
more people. 

And Julian.

## License
GPL 3.0? 

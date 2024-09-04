# TuMag's pipeline

Tools for reading and processing raw data from TuMag.

Brief description: 
- Check_image_id.py : Python function to be called throuth termiinal to read images through ID.
- config.py : Config fie with headers info, observation modes info, etc
- demodulation.py : Module to compute demodulation. 
- field_stop_finder.py : Module to find field stop and align images. 
- fringes.py : Module to clean images from fringes. 
- image_handler.py : Module to read .img files, process observation modes, flat-field modes.
- master_dark.py : Module to compute dark current. 
- master_flat_field : Module to compute master flat-fields- 
- organizer.py : Python script to generate the IDs for the images.   
- requirements.txt : Dependencies to run the pipeline.
- utils.py : Helper module to process raw images. 
- vlos.py : Module to compute the cog velocities. 
- xtalk.py : Module to compute and correct x-talk. 

Disclaimer: Very early still...

## Useful guides
- [Installation](Documents/Installation.md)
- [Working with TuMag's data](Documents/Working_with_Tumags_data.md)

## Overview of TuMag's observations

### Timelines

TuMag observations were stuctured through  **timelines**, a list of pre-configured commands where calibration and scientific observations were already programmed. 

A summary of the observations and the corresponding timelines can be foound in [Observation summary](Documents/Sunrise_summary_observations.pdf). 

For a guide describing TuMag's commands for each timeline see [Julián's webpage](https://www.uv.es/jublanro/TuMag_timeline_reference.html). 

<ins>Important:</ins> When working with the data have in mind that the timelines were paused and stopped, some commands were run manually, others aborted, etc ... The webpage should be used as a reference. 

All observations (scientific and calibrations) are recorded and grouped in their corresponding Timeline in the [Tumag's Logbook](Documents/TuMagCompass.csv). 

Tumag observations are operated through pre-configured **observation modes**. 

| **Observing Mode** | **Filter** | **N_{$\lambda$}** | **Modulation scheme** | 
|----------|----------|----------|----------|
| 0p    | 517 nm   | 15  | Vectorial   | 
| 0s    | 517 mn   | 15  | Vectorial (1 Mod)   |
| 1    | 517 nm   | 10 | Vectorial  | 
| 2.02    | 525.02 nm  | 8  | Vectorial  |
| 2.06    | 525.06  | 8  | Vectorial  |
| 3.02    | 525.02 nm  | 5  | Longitudinal  |
| 3.06    | 525.06  | 5  | Longitudinal  |
| 4    | 517  | 3  | Data 19  | Vectorial  |
| 5.02    | 525.02 nm  | 3  | Vectorial  |
| 5.06    | 525.06  | 3  | Vectorial  |

## Contributing

SPG team.

## License
GPL 3.0? 

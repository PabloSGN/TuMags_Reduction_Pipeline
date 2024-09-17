# Module guide

Brief description: 
- Check_image_id.py : Python function to be called through terminal to read images through ID.
- config.py : Config file with headers info, observation modes info, etc
- demodulation.py : Module to compute demodulation. 
- field_stop_finder.py : Module to find field stop and align images. 
- fringes.py : Module to clean images from fringes. 
- image_handler.py : Module to read .img files, process observation modes, flat-field modes. [Documentation](./Modules_documentations/image_handler.md)
- master_dark.py : Module to compute dark current. [Documentation](./Modules_documentations/master_dark.md)
- master_flat_field : Module to compute master flat-fields. [Documentation](./Modules_documentations/master_flatfield.md)
- organizer.py : Python script to generate the IDs for the images.   
- requirements.txt : Dependencies to run the pipeline.
- utils.py : Helper module to process raw images. 
- vlos.py : Module to compute the cog velocities. 
- xtalk.py : Module to compute and correct x-talk. 
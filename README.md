# TuMag's pipeline

Tools for reading and processing raw data from TuMag.

Disclaimer: Very early still...

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Installation of required libraries. 

### IAA's server guide

Install the required libraries through a virtual enviroment:

- Anywhere on your directories -> run the command: "python3 -m venv $myenv$", substituting $myenv$ for whatever name you want to give it. 
- Activate the enviroment: run the command: "source $myenv$/bin/activate.csh", again substitute the name. 
- Install dependencies:  "pip install -r requirements.txt"

You should be able to run any script after the installation finishes with the enviroment activated. 

## Usage

### Working with image ID's

This option is only valid when working with IAA's server. 

Image IDs have the format DXX-YYYY. The DXX points to the csv file in the Organized files directory where the image's path is stored. And the YYYY is the individual index of the image within that file (first column). 
Image IDs can be found in the csv files. A guide to the observations labeled with image indexes can be found in the Documents directory under the name of: "TuMag's Compass". 

- Reading a specific image. 

The "image_handler" module incñudes the function: "read_ID". Provided with an index of an image returns the Image (2D numpy array) and its header. 
Example of use: 
    I, H = image_handler.read_ID("D11-1917")
Please refer to the function for more info. 

- Obtaining the paths of multiple images

The "image_handler" module includes the function: "get_images_paths". This function returns the paths of the images between two indexes of a specific file. It can be given multiple queries if images from different files are required.
Example of use: 
    paths = image_handler([["D11", 1917, 5620], ["D12", 302, 17000]])

    The paths variable will include all paths for the images between index 1917 and 5620 of file "D11.csv" and from 302 to 17000 of file "D12.csv"


## Contributing
- Pablo (psanta@iaa.es)


## License
GPL 3.0? 
# TuMag's pipeline

Tools for reading and processing raw data from TuMag.

Disclaimer: Very early still...

## Table of Contents
- [Installation](#installation-of-required-libraries)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Installation of required libraries. 

All the required libraries are included in the [requirements.txt](requirements.txt) file. Install them globally or through an enviroment. 

If using the IAA server (recomended to be able to use images IDs) the use of an enviroment is mandatory. 

### Installing the dependencies with an enviroment. 

- Anywhere on your directories -> run the command: 
```shell
"python3 -m venv $myenv$"
```
Substituting "$myenv$" for whatever name you want to give it. 
- Activate the enviroment: run the command: 
```shell
source $myenv$/bin/activate.csh
```
Again substitute the name of "$myenv$". 
- Install dependencies: 
```shell
pip install -r requirements.txt
```

You should be able to run any script after the installation finishes with the enviroment activated. 

## Usage

### Working with image ID's

This option is only valid when working with IAA's server. 

Image IDs have the format DXX-YYYY. The DXX points to the csv file in the [Organized files](Organized_files/) directory where the image's path is stored. And the YYYY is the individual index of the image within that file (first column). 
Image IDs can be found in the csv files. A guide to the observations labeled with image indexes can be found in the Documents directory under the name of: "TuMag's Compass" (to be uploaded). 

#### <ins>Reading a specific image.</ins>

The [image_handler](image_handler.py) module includes the function: [read_ID](image_handler.py#L258). Provided with an index of an image returns the Image (2D numpy array) and its header. 

- Example of use:
```python
I, H = image_handler.read_ID("D11-1917")
```
Please refer to the function for more info. 

#### <ins>Obtaining the paths of multiple images.</ins>

The [image_handler](image_handler.py) module includes the function: [get_images_paths](image_handler.py#L216). This function returns the paths of the images between two indexes of a specific file. It can be given multiple queries if images from different files are required.

- Example of use:
```python
paths = image_handler.get_images_paths([["D11", 1917, 5620], ["D12", 302, 17000]])
```

The paths variable will include all paths for the images between index 1917 and 5620 of file "D11.csv" and from 302 to 17000 of file "D12.csv"

#### <ins>Reading a specific image from terminal:</ins>

The module [Check_image_id](Check_image_id.py) can be called in the terminal line to read an image given an image index. By default it only prints the observation mode and counter. But if the word "plot" and/or "header" is added to the command line will also plot and print the header, respectively.

- Example of use:
```shell
python3 Check_image_id.py D11-1500 plot header
```


## Contributing
- Pablo (psanta@iaa.es)


## License
GPL 3.0? 
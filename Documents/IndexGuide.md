# Working with TuMag's data. 

TuMag data is stored in IAA's server in "work/obs/TuMAG_data" and stored in different directories depending on the the day and hour of the observation. 

To identify the images belonging to scientific observations, all observation-related images of TuMag have been assigned a unique ID. The ID's have the structure: DXX-YYYY. Where the DXX stands for the day of the measurement: (D09, D10, D11, ... , D16), and YYYY is a cardinal number to identify the specific image within that day. 

There is a csv file for each day where all images are labeled with the specific path to the image within the server. The files are provided in this repository in the the [Organized files](../Organized_files/) folder. 

To identify an image just go to the csv file corresponding to the "DXX" and search the image labeled with YYYY. 

#### <ins>Example</ins>

Image: D11-102. 

From file [D11.csv](../Organized_files/D11.csv#L103) in line 103 (ID's start at 0):
- 102,/work/obs/TuMAG_data/2024_07_11_00_00_00_123/2024_07_11_00_00_00_123/2024_07_11_00_00_45_321_1_9189.img


## Identifying observations. 

All observations are recorded in the [TuMag's Logbook](TuMagCompass.csv) file. All timelines are detailed in the file with all the observation and calibration blocks indicated. 

Every observation block has a column where the ID's are noted down. The format of labelling the observation block is "DXX-YYYY-ZZZZ". This indicates that the images composing the observation range from the image DXX-YYYY to DXX-ZZZZ, including all images in between. Some times different quearys appear in these cells, meaning that the images are separated.

Please read the remarks (sometimes difficult to understand even for the one who wrote them...)because sometimes are missing files, aborted observations, etc...   

#### <ins>Example of labeled timeline</ins>
![Screenshot](Images_for_guides/example_guide_indexes.png)
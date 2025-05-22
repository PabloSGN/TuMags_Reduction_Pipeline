# ---------------------------- DESCRIPTION --------------------------------------- #
"""

Module to read an image from terminal using the image ID (IAA's server) 

run command:

python/python3 Check_image_id.py image_ID 

Some key words can be added to perform additional tasks:
     - plot -> Plots the image
    - header -> Reads the header.
    - bin -> Bins the image for  faster reading.

Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

import sys
import image_handler as ih

if __name__ == "__main__":

    args = sys.argv

    if "plot" in args:
        plotflag = True
    else:
        plotflag = False

    if "header" in args:
        headerflag = True
    else:
        headerflag = False

    if "bin" in args:
        binningflag = True
    else:
        binningflag = False


    I, H = ih.read_ID(args[1], verbose = True, plotflag=plotflag, 
                      header = headerflag, binning = binningflag)
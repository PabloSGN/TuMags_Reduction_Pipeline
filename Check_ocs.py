# ---------------------------- DESCRIPTION --------------------------------------- #
"""

Module to read a series of images and separate them with respect to their OC.

run command:

python/python3 Check_ocs.py image_IDs 

image_IDs Format -> DXX-init-end

Some key words can be added to perform additional tasks:
    - flat -> Uses the flat mode (expects more images)
    - save -> Saves the result in a pickle file of name (OCs.pickle) for later use.

Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

import sys
import os
import pickle
import image_handler as ih

if __name__ == "__main__":

    args = sys.argv

    if "flat" in args:
        flatflag = True
    else:
        flatflag = False
   
    paths = ih.get_images_paths(args[1])
    ocs = ih.separate_ocs(paths, verbose = True, flat_fieldmode = flatflag)

    if "save" in args:
        count = 1

        filename = "OCs.pickle"
        while os.path.exists(filename):
            filename = f"OCs_{count}.pickle"
            count += 1
        
        with open(filename, "wb") as output_file:
            pickle.dump(ocs, output_file)

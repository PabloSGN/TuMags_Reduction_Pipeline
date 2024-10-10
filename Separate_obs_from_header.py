import sys
import os
import pickle

import image_handler as ih

if __name__ == "__main__":

    args = sys.argv

    ocs = ih.obs_mode_separator(ih.get_images_paths(args[1]), verbose = True)
    
    if "save" in args:
        count = 1

        filename = "OCs_headers.pickle"
        while os.path.exists(filename):
            filename = f"OCs_headers_{count}.pickle"
            count += 1
        
        with open(filename, "wb") as output_file:
            pickle.dump(ocs, output_file)
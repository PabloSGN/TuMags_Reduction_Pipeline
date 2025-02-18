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

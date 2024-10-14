import sys
import os
import pickle
import image_handler as ih

if __name__ == "__main__":

    args = sys.argv
   
    paths = ih.get_images_paths(args[1])

    paths = [x for x in paths if "_0_" in x]

    for x in paths[5725:5730]:
        print(os.path.basename(x))
   
    print(paths[0])
    ih.check_timestamps(paths, verbose = True)



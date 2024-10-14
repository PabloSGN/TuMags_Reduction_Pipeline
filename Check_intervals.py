import sys
import os
import pickle
import image_handler as ih

if __name__ == "__main__":

    args = sys.argv
   
    paths = ih.get_images_paths(args[1])
    ih.check_timestamps(paths, verbose = True)



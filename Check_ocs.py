import sys
import image_handler as ih

if __name__ == "__main__":

    args = sys.argv
   
    paths = ih.get_images_paths(args[1])
    ocs = ih.separate_ocs(paths, verbose = True)
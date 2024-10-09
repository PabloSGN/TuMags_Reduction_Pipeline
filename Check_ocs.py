import sys
import image_handler as ih

if __name__ == "__main__":

    args = sys.argv

    if "flat" in args:
        flatflag = True
    else:
        flatflag = False
   
    paths = ih.get_images_paths(args[1])
    ocs = ih.separate_ocs(paths, verbose = True, flat_fieldmode = flatflag)
import sys
import image_handler as ih

if __name__ == "__main__":

    args = sys.argv
    I, H = ih.read_ID(args[0], int(args[1]), verbose = True)
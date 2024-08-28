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

    I, H = ih.read_ID(args[1], verbose = True, plotflag=plotflag, header = headerflag)
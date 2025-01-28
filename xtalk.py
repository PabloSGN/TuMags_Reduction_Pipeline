# ---------------------------- DESCRIPTION --------------------------------------- #

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in 
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mp
from astropy.io import fits
from scipy.optimize import curve_fit

# ------------------------------ CONFIG ------------------------------------------ #

plt.style.use("default")
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"

# ------------------------------  CODE  ------------------------------------------ # 

def compute_xtalk_full(data, center_flag = False, savefig = False, figname = "xtalk.png", rad = 250):

        xtalked = np.copy(data)

        fig, axs = plt.subplots(1, 3, figsize = (15, 5))

        xtalked[1], axq = compute_xtalk(data, "I", "Q", fig = fig, ax = axs[0], center = center_flag, rad = rad)
        xtalked[2], axu = compute_xtalk(data, "I", "U", fig = fig, ax = axs[1], center = center_flag, rad = rad)
        xtalked[3], axv = compute_xtalk(data, "I", "V", fig = fig, ax = axs[2], center = center_flag, rad = rad)
        plt.tight_layout()
        if savefig:
            plt.tight_layout()
            plt.savefig(figname, bbox_inches = "tight")


        return xtalked

def compute_xtalk(input, stokes_x, stokes_y, fig, ax, center = False, rad = 250):
        
        ind = {"I" : 0,
               "Q" : 1,
               "U" : 2,
               "V" : 3}
        
        ndim = np.shape(input)[-1]

        if center:
                cent_coords = ndim // 2
                data = input[:, cent_coords - rad : cent_coords + rad, 
                                cent_coords - rad : cent_coords + rad]               
        else:
              data = np.copy(input)

        # Filtering out of the fit bright and dark points
        I = data[0]
        norm = np.median(I)
        points = (I > 0.7 * norm) & (I < 1.3 * norm)

        #Selecting points for the fit. 
        xaxis = data[ind[stokes_x]][points].flatten() #/ norm
        yaxis = data[ind[stokes_y]][points].flatten() #/ norm

        def line(x, a, b):
            return x * a + b

        # Fitting
        fit_params = curve_fit(line, xaxis, yaxis)
        coeffs = fit_params[0]
        fit = np.poly1d(coeffs)

        # Generate the 2D histogram for representation
        ax.hist2d(xaxis, yaxis, bins=(50, 50), cmap=plt.cm.jet)

        ax.set_title(f"{stokes_x} -> {stokes_y}")
        ax.plot([np.min(xaxis), np.max(xaxis)], [fit(np.min(xaxis)), fit(np.max(xaxis))], lw = 3, ls = '-', c = 'w', 
                label = f"Fit (a: {round(coeffs[0], 3)})")

        ax.set_xlabel(stokes_x)
        ax.set_ylabel(stokes_y)

        ax.legend(edgecolor = 'k', facecolor = 'w')

        corr = input[ind[stokes_y]] - (input[ind[stokes_x]] * coeffs[0] - coeffs[1]) 

        return corr, ax
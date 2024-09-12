# ---------------------------- DESCRIPTION --------------------------------------- #

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in 
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mp
from astropy.io import fits
from scipy.optimize import curve_fit
import mpl_scatter_density
from scipy.stats import binned_statistic


# ------------------------------ CONFIG ------------------------------------------ #

plt.style.use("default")
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"


# Color map for scatter plot representation
white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
        (0, '#000000'),
        (1e-20, '#330033'),
        (0.2, '#990066'),
        (0.4, '#CE22CE'),
        (0.6, '#CC3333'),
        (0.8, '#FF6600'),
        (1, '#FFF700'),
        ], N=256)

# ------------------------------  CODE  ------------------------------------------ # 

def compute_xtalk(data, stokes_x, stokes_y, savefig = False, filename = 'image'):
        
        ind = {"I" : 0,
               "Q" : 1,
               "U" : 2,
               "V" : 3}
        
        
        I = data[0]
        norm = np.median(I[300:-300, 300:-300])
        points = (I > 0.8 * norm) & (I < 1.2 * norm)

        norm = np.median(I[300:-300, 300:-300])

        xaxis = data[ind[stokes_x]][points].flatten() #/ norm
        yaxis = data[ind[stokes_y]][points].flatten() #/ norm

        print(np.min(xaxis))

        def line(x, a, b):
            return b + x * a

        fit_params = curve_fit(line, xaxis, yaxis )

        coeffs = fit_params[0]


        print(coeffs)

        fit = [np.min(xaxis) * coeffs[0] + coeffs[1], np.max(xaxis) * coeffs[0] + coeffs[1]]
        
        
        def using_mpl_scatter_density(fig, x, y):
            ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
            density = ax.scatter_density(x, y, cmap=white_viridis)
            fig.colorbar(density, label='Number of points per pixel')
            
            return ax
        

        line = [np.min(xaxis) * fit[0] + fit[1], np.max(xaxis) * fit[0] + fit[1]]

        fig = plt.figure(figsize = (5.5, 5))
        ax = using_mpl_scatter_density(fig, xaxis, yaxis)


        ax.set_title(f"{stokes_x} -> {stokes_y}")
        stats = binned_statistic(xaxis, yaxis, statistic = 'mean', bins = 50)
        #ax.plot(xaxis, stats[0], lw = 3, ls = '-', c = 'dodgerblue', label = 'Mean value')
        ax.plot([np.min(xaxis), np.max(xaxis)], fit, lw = 3, ls = '-', c = 'w', 
                label = f"Fitting to line: f(x) = ax + b (a: {round(coeffs[0], 3)})")

        #ax.set_xlim(np.min(xaxis), np.max(xaxis))

        ax.set_xlabel(stokes_x)
        ax.set_ylabel(stokes_y)

        ax.legend(edgecolor = 'k', facecolor = 'darkorange')

        corr = data[ind[stokes_y]] - (coeffs[1] + data[ind[stokes_x]] * coeffs[0])      
        return corr
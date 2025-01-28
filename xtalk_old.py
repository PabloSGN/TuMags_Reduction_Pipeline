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
def reshape_into_16_quadrants(data):
    # Reshape the last two dimensions into a 4x4 grid of (354, 354)
    reshaped = data.reshape(4, 4, 354, 4, 354)
    # Rearrange axes to group quadrants into a single dimension
    return reshaped.transpose(0, 1, 3, 2, 4).reshape(4, 16, 354, 354)

def reverse_quadrants(quadrants):
# Reverse the reshaping process
    reshaped = quadrants.reshape(4, 4, 4, 354, 354).transpose(0, 1, 3, 2, 4)
    reconstructed = reshaped.reshape(4, 1416, 1416)

    return reconstructed


def compute_xtalk_quadrants(data, stokes_x, stokes_y, fig, fig_ind):
      
        ind = {"I" : 0,
                "Q" : 1,
                "U" : 2,
                "V" : 3}

        I = data[0]
        norm = np.median(I)
        points = (I > 0.5 * norm) & (I < 1.5 * norm)

        xaxis = data[ind[stokes_x]][points].flatten() #/ norm
        yaxis = data[ind[stokes_y]][points].flatten() #/ norm

        print(np.min(xaxis))

        def line(x, a, b):
                return b + x * a

        fit_params = curve_fit(line, xaxis, yaxis )

        coeffs = fit_params[0]

        print(coeffs)

        fit = [np.min(xaxis) * coeffs[0] + coeffs[1], np.max(xaxis) * coeffs[0] + coeffs[1]]

        def using_mpl_scatter_density(fig, fig_ind,  x, y):
                ax = fig.add_subplot(4, 4, fig_ind + 1, projection='scatter_density')
                ax.axis("off")
                density = ax.scatter_density(x, y, cmap=white_viridis)
                #fig.colorbar(density, label='Number of points per pixel')
                
                return ax

        ax = using_mpl_scatter_density(fig, fig_ind, xaxis, yaxis)

        line = [np.min(xaxis) * fit[0] + fit[1], np.max(xaxis) * fit[0] + fit[1]]

        stats = binned_statistic(xaxis, yaxis, statistic = 'mean', bins = 50)
        #ax.plot(xaxis, stats[0], lw = 3, ls = '-', c = 'dodgerblue', label = 'Mean value')
        ax.plot([np.min(xaxis), np.max(xaxis)], fit, lw = 3, ls = '-', c = 'w', 
                label = f"Fitting to line: f(x) = ax + b (a: {round(coeffs[0], 3)})")

        print(data.shape)
        print(ind)
        print(stokes_y, stokes_x)
        print(ind[stokes_y], ind[stokes_x])

        corr = data[ind[stokes_y]] - (coeffs[1] + data[ind[stokes_x]] * coeffs[0]) 

        return corr, ax
      

def compute_xtalk_full_quadrants(data):
        
        quads = reshape_into_16_quadrants(data)
        xtalked_quads = np.copy(quads)

        figq = plt.figure(figsize = (0, 10))
        figu = plt.figure(figsize = (0, 10))
        figv = plt.figure(figsize = (0, 10))

        for qq in range(16):
                xtalked_quads[1, qq], axq = compute_xtalk_quadrants(quads[:, qq], "I", "Q", fig = figq, fig_ind = qq)
                xtalked_quads[2, qq], axu = compute_xtalk_quadrants(quads[:, qq], "I", "U", fig = figu, fig_ind = qq)
                xtalked_quads[3, qq], axv = compute_xtalk_quadrants(quads[:, qq], "I", "V", fig = figv, fig_ind = qq)

        figq.suptitle("I -> Q")
        figu.suptitle("I -> U")
        figv.suptitle("I -> V")

        """figq.tight_layout()
        figu.tight_layout()
        figv.tight_layout()"""

        return reverse_quadrants(xtalked_quads)

def compute_xtalk_full(data, center = False):

        xtalked = np.copy(data)

        fig = plt.figure(figsize = (19, 5))

        xtalked[1], axq = compute_xtalk(data, "I", "Q", fig = fig, fig_ind = 1, center = center)
        xtalked[2], axu = compute_xtalk(data, "I", "U", fig = fig, fig_ind = 2, center = center)
        xtalked[3], axv = compute_xtalk(data, "I", "V", fig = fig, fig_ind = 3, center = center)

        return xtalked


def compute_xtalk(input, stokes_x, stokes_y, fig, fig_ind, center = False):
        
        ind = {"I" : 0,
               "Q" : 1,
               "U" : 2,
               "V" : 3}
        
        ndim = np.shape(input)[-1]

        if center:
                cent_coords = ndim // 2
                print(cent_coords)
                data = input[:, cent_coords - 100 : cent_coords + 100, cent_coords - 100 : cent_coords + 100]               
        else:
              data = np.copy(input)
              

        I = data[0]
        norm = np.median(I)
        points = (I > 0.7 * norm) & (I < 1.3 * norm)

        xaxis = data[ind[stokes_x]][points].flatten() #/ norm
        yaxis = data[ind[stokes_y]][points].flatten() #/ norm

        print(np.min(xaxis))

        def line(x, a, b):
            return b + x * a

        fit_params = curve_fit(line, xaxis, yaxis )

        coeffs = fit_params[0]

        print(coeffs)

        fit = [np.min(xaxis) * coeffs[0] + coeffs[1], np.max(xaxis) * coeffs[0] + coeffs[1]]
        
        def using_mpl_scatter_density(fig, fig_ind,  x, y):
            ax = fig.add_subplot(1, 3, fig_ind, projection='scatter_density')
            density = ax.scatter_density(x, y, cmap=white_viridis)
            fig.colorbar(density, label='Number of points per pixel')
            
            return ax

        ax = using_mpl_scatter_density(fig, fig_ind, xaxis, yaxis)
        
        line = [np.min(xaxis) * fit[0] + fit[1], np.max(xaxis) * fit[0] + fit[1]]

        ax.set_title(f"{stokes_x} -> {stokes_y}")
        stats = binned_statistic(xaxis, yaxis, statistic = 'mean', bins = 50)
        #ax.plot(xaxis, stats[0], lw = 3, ls = '-', c = 'dodgerblue', label = 'Mean value')
        ax.plot([np.min(xaxis), np.max(xaxis)], fit, lw = 3, ls = '-', c = 'w', 
                label = f"Fitting to line: f(x) = ax + b (a: {round(coeffs[0], 3)})")

        #ax.set_xlim(np.min(xaxis), np.max(xaxis))

        ax.set_xlabel(stokes_x)
        ax.set_ylabel(stokes_y)

        ax.legend(edgecolor = 'k', facecolor = 'darkorange')

        corr = input[ind[stokes_y]] - (coeffs[1] + input[ind[stokes_x]] * coeffs[0]) 

        return corr, ax
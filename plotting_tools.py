# ---------------------------- DESCRIPTION --------------------------------------- #

"""
Plotting tools
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in libs
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap

# Own-libs
import alignment as al
# ------------------------------ CONFIG ------------------------------------------ #
def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)
# --------------------------------------------------------------------------------- # 

def plot_mosaic(data, clims = None, norm = None):

    # Data shape -> 4 (mods / stokes ) x quads x Nx x Ny

    fig, axs = plt.subplots(8, 8, figsize = (16.5, 15), gridspec_kw={'wspace': 0, 'hspace': 0})

    if norm is None:
        norm = np.mean(data)

    flat_ax = axs.flatten()

    for i, ax in enumerate(flat_ax):
        ax.set_xticks([])
        ax.set_yticks([])

    q1 = np.array([0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27])
    q2 = q1 + 4
    q3 = q1 + 32
    q4 = q1 + 36

    quad_coords = [q1, q2, q3, q4]

    
    flag1 = True
    flag2 = True
    flag3 = True
    flag4 = True
    for i in range(4):
        coords = quad_coords[i]
        
        for quad_ind, ax_ind in enumerate(quad_coords[i]):
            
            if quad_ind == 0 and flag1:
                if clims is None:
                    im1 = flat_ax[ax_ind].imshow(data[i, quad_ind] / norm, clim = (0.5, 1.5), cmap = "gray")
                else:
                    im1 = flat_ax[ax_ind].imshow(data[i, quad_ind] / norm, clim = clims[i], cmap = "gray")
                flag1 = False
            elif quad_ind == 0 and flag2:
                if clims is None:
                    im2 = flat_ax[ax_ind].imshow(data[i, quad_ind] / norm, clim = (0.5, 1.5), cmap = "gray")
                else:
                    im2 = flat_ax[ax_ind].imshow(data[i, quad_ind] / norm, clim = clims[i], cmap = "gray")
                flag2 = False
            elif quad_ind == 0 and flag3:
                if clims is None:
                    im3 = flat_ax[ax_ind].imshow(data[i, quad_ind] / norm, clim = (0.5, 1.5), cmap = "gray")
                else:
                    im3 = flat_ax[ax_ind].imshow(data[i, quad_ind] / norm, clim = clims[i], cmap = "gray")
                flag3 = False

            elif quad_ind == 0 and flag4:
                if clims is None:
                    im4 = flat_ax[ax_ind].imshow(data[i, quad_ind] / norm, clim = (0.5, 1.5), cmap = "gray")
                else:
                    im4 = flat_ax[ax_ind].imshow(data[i, quad_ind] / norm, clim = clims[i], cmap = "gray")
                flag4 = False
            else:
                if clims is None:
                    flat_ax[ax_ind].imshow(data[i, quad_ind] / norm, clim = (0.5, 1.5), cmap = "gray")
                else:
                    flat_ax[ax_ind].imshow(data[i, quad_ind] / norm, clim = clims[i], cmap = "gray")

    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    # Add the colorbars
    # Define positions: (x, y, width, height) as fractions of the figure
    left_cb1_pos = [0.02, 0.55, 0.02, 0.4]  # Top-left
    left_cb2_pos = [0.02, 0.05, 0.02, 0.4]   # Bottom-left
    right_cb1_pos = [0.95, 0.55, 0.02, 0.4]  # Top-right
    right_cb2_pos = [0.95, 0.05, 0.02, 0.4]   # Bottom-right

    # Create colorbars
    cbar1 = fig.colorbar(im1, cax=fig.add_axes(left_cb1_pos), orientation='vertical')
    cbar2 = fig.colorbar(im3, cax=fig.add_axes(left_cb2_pos), orientation='vertical')
    cbar3 = fig.colorbar(im2, cax=fig.add_axes(right_cb1_pos), orientation='vertical')
    cbar4 = fig.colorbar(im4, cax=fig.add_axes(right_cb2_pos), orientation='vertical')

    # Set ticks and labels for left colorbars
    for cbar in [cbar1, cbar2]:
        cbar.ax.yaxis.set_ticks_position('left')  # Ticks on the left
        cbar.ax.yaxis.set_label_position('left')  # Label on the left


    plt.tight_layout(rect=[0.04, 0, 0.96, 1])  # Adjust layout to fit colorbars

def plot_quadrant(quadrant, norm = None, clims = None):
    
    fig, axs = plt.subplots(2, 2, figsize = (12, 10))

    if norm is None:
        norm = 1

    for ind, ax in enumerate(axs.flatten()):

        if clims is not None:
            im = ax.imshow(quadrant[ind] / norm, cmap = "gray", clim = clims[ind])
        else:
            im = ax.imshow(quadrant[ind] / norm, cmap = "gray")
        colorbar(im)

        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
import os  
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt 

import matplotlib.animation as manimation
import numpy as np
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Hi David, 

# thank you for sharing this. I have downloaded and accessed it. 
# The provided package looks already quite good. However, there are some details to tackle:

# The reference images look nice, however, it would be great if you could select the two most representative ones. I would suggest a continuum map from the beginning and a selected highlight (e.g. a Stokes V map, if available) from one line. I'd suggest <OBS_ID>_continuum.png and <OBS_ID>_highlight.png as naming template. With this I can just automatically insert the files to the website for all observations.
# The plots should include as little white space as possible.

# As you can see on the observation-overview (https://sr3test.mps.mpg.de/obs) page, the thumbnail isn't all that big. However, it seems a good idea to include the observed line in the plot title... 
# What do you think of the following: Change the title of the plots to include the observed line. Additionally, create a thumbnail version of the <OBS_ID>_continuum.png file, without any title, axis labels or color bar.
# The thumbnail will then be displayed in the list overview, while the two others will be used in the observation details view (https://sr3test.mps.mpg.de/obs/01_QSUN) on the left. 

# For the details plots:
# It would be best if this would be one panel / one move. This would allow the user to play all panels simultaneously and give you the control over their interaction. 
# Please note that landscape or wide-screen aspect is best for the place reserved here. So for 01_QSUN it seems best to have a four-axis plot with Fe (I), Fe (V), Mg (I), MG (V) in a row.
# Also, it allows to have one PNG and one MOVIE for each observation. As some obs might only have one or two panels. To cope with that, the tool would need to know how many plots for what obs. That is too much detail knowledge to implement. This should be called <OBS_ID>_<TUMAG_ID>_details.(png/mp4) (or similar, as long as it allows for automatic connection to a tumag-obs.)

# For the fits I'll ask Andi  to have a look if his data download request processor can cope with these files and how to streamline the format between the instruments. 
# I also set Francisco in CC, if a discussion on the fits format between instruments is needed. 
# @Andi: Can you have a look at the fits provided? 

# Cheers,
# Johannes

name = '01_QSUN'
line =  'Fe2.02'#'Fe1' # "Fe2.02"
clim_V = (-1, 1)
clim_I = (0.6, 1.25)
roi = [0, -1, 0, -1]
wave_I = 0
wave_V = 2

line =  'Mg' #'Fe2.02'#'Fe1' # "Fe2.02"
clim_V = (-1, 1)
clim_I = (0.6, 1.25)
roi = [0, -1, 0, -1]
wave_I = -1
wave_V = 2

vers = 'v0.2'# '1.1'
dir = './' #'/Volumes/TuMag_d1/reduccion/'+name

files = sorted(os.listdir(dir))
files = [f for f in files if f.endswith('.fits') and not f.startswith('._')]
files = [f for f in files if vers in f ]
files = [f for f in files if line in f ]
# files = [f for f in files if "Fe2.02" in f ]
filed = files.sort()
files = [os.path.join(dir, f) for f in files]
if line == 'Mg':
    files = files[1:]
print(files)
for f in files:
    header = fits.getheader(f)

# -------------------------
# Parámetros
# -------------------------

output_dir = "reference_pngs"
os.makedirs(output_dir, exist_ok=True)




# ---------------------------------
# Función auxiliar para plottear
# ---------------------------------
def plot_stokes(image, clim, title, cbar_label, outname, date_obs):
    fig, ax = plt.subplots(figsize=(4, 4))

    im = ax.imshow(image, cmap="gray", clim=clim)
    ax.set_title(f"{title}\n{date_obs}", fontsize=10)
    ax.axis("off")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    cbar = plt.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(cbar_label, fontsize=9)

    plt.tight_layout()
    plt.savefig(outname, dpi=150)
    plt.close()

# =================================
# PRIMER FRAME
# =================================
hdu = fits.open(files[0])
data = hdu[0].data
header = hdu[0].header
date_obs = header.get("DATE-OBS", "UNKNOWN")

I = data[wave_I, 0, roi[0]:roi[1], roi[2]:roi[3]]
V = data[wave_V, 3, roi[0]:roi[1], roi[2]:roi[3]] * 100.0
print(data.shape)
plot_stokes(
    I, clim_I,
    "Stokes I",
    "[Intensity]",
    f"{output_dir}/first_I_{name}_{line}.png",
    date_obs
)

plot_stokes(
    V, clim_V,
    "Stokes V",
    "[%]",
    f"{output_dir}/first_V_{name}_{line}.png",
    date_obs
)

hdu.close()

# =================================
# ÚLTIMO FRAME
# =================================
hdu = fits.open(files[-1])
data = hdu[0].data
header = hdu[0].header
date_obs = header.get("DATE-OBS", "UNKNOWN")

I = data[wave_I, 0, roi[0]:roi[1], roi[2]:roi[3]]
V = data[wave_V, 3, roi[0]:roi[1], roi[2]:roi[3]] * 100.0

plot_stokes(
    I, clim_I,
    "Stokes I",
    "[Intensity]",
    f"{output_dir}/last_I_{name}_{line}.png",
    date_obs
)

plot_stokes(
    V, clim_V,
    "Stokes V",
    "[%]",
    f"{output_dir}/last_V_{name}_{line}.png",
    date_obs
)

hdu.close()


# MOVIE

# Directorio de salida
out_dir = f"movies"
os.makedirs(out_dir, exist_ok=True)

# -------------------------
# Loop sobre los FITS

# -------------------------
# Directorios de salida
# -------------------------
out_dir_I = f"png_I_{name}_{line}"
out_dir_V = f"png_V_{name}_{line}"

os.makedirs(out_dir_I, exist_ok=True)
os.makedirs(out_dir_V, exist_ok=True)

# -------------------------
# Loop sobre los FITS
# -------------------------
for idx, f in enumerate(files):
    print(f"Frame {idx+1}/{len(files)}")

    # Datos y cabecera
    hdu = fits.open(f)
    data = hdu[0].data
    header = hdu[0].header

    date_obs = header.get("DATE-OBS", "UNKNOWN")

    # ---------
    # Stokes I
    # ---------
    I = data[wave_I, 0, roi[0]:roi[1], roi[2]:roi[3]]

    fig, ax = plt.subplots(figsize=(4, 4))
    im = ax.imshow(I, cmap="gray", clim=clim_I)
    ax.set_title(f"Stokes I\n{date_obs}", fontsize=10)
    ax.axis("off")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('[Ic]', fontsize=9)   # ajusta unidades si quieres

    plt.tight_layout()

    fname_I = (
        f"{out_dir_I}/"
        f"{name}_{line}_I_wave{wave_I}_"
        f"{date_obs.replace(':','-')}_"
        f"frame_{idx:04d}.png"
    )

    plt.savefig(fname_I, dpi=150)
    plt.close()

    # ---------
    # Stokes V
    # ---------
    V = data[wave_V, 3, roi[0]:roi[1], roi[2]:roi[3]]

    fig, ax = plt.subplots(figsize=(4, 4))
    im = ax.imshow(V*100, cmap="gray", clim=clim_V)
    ax.set_title(f"Stokes V\n{date_obs}", fontsize=10)
    ax.axis("off")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('[%]', fontsize=9)

    plt.tight_layout()

    fname_V = (
        f"{out_dir_V}/"
        f"{name}_{line}_V_wave{wave_V}_"
        f"{date_obs.replace(':','-')}_"
        f"frame_{idx:04d}.png"
    )

    plt.savefig(fname_V, dpi=150)
    plt.close()

    hdu.close()


fps = 2

script_dir = os.path.dirname(os.path.abspath(__file__))
ffmpeg_path = os.path.join(script_dir, 'ffmpeg')

plt.rcParams['animation.ffmpeg_path'] = ffmpeg_path
FFMpegWriter = manimation.writers['ffmpeg']

writer = FFMpegWriter(fps=fps)

# -------------------------------------------------
# Figura inicial (UNA sola vez)
# -------------------------------------------------
fig, ax = plt.subplots(figsize=(4, 4))

# Imagen "dummy" para inicializar la colorbar
dummy = ax.imshow(
    fits.getdata(files[0])[wave_I, 0, roi[0]:roi[1], roi[2]:roi[3]],
    cmap="gray",
    clim=clim_I
)

ax.axis("off")

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

cbar = plt.colorbar(dummy, cax=cax)
cbar.ax.tick_params(labelsize=8)

# -------------------------------------------------
# Animación
# -------------------------------------------------

with writer.saving(fig, f"{out_dir}/movie_{name}_{line}_StokesI.mp4", dpi=120):
    for f in files:
        hdu = fits.open(f)
        data = hdu[0].data
        header = hdu[0].header
        date_obs = header.get("DATE-OBS", "UNKNOWN")

        I = data[wave_I, 0, roi[0]:roi[1], roi[2]:roi[3]]

        dummy.set_data(I)
        ax.set_title(f"Stokes I {date_obs}", fontsize=10)

        writer.grab_frame()
        hdu.close()

plt.close(fig)

# -------------------------------------------------
# Figura inicial (UNA sola vez)
# -------------------------------------------------
fig, ax = plt.subplots(figsize=(4, 4))

# Imagen "dummy" para inicializar la colorbar
dummy = ax.imshow(
    fits.getdata(files[0])[wave_V, 3, roi[0]:roi[1], roi[2]:roi[3]],
    cmap="gray",
    clim=clim_V
)

ax.axis("off")

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

cbar = plt.colorbar(dummy, cax=cax)
cbar.ax.tick_params(labelsize=8)

# -------------------------------------------------
# Animación
# -------------------------------------------------
fig.subplots_adjust(right=0.85)
with writer.saving(fig, f"{out_dir}/movie_{name}_{line}_StokesV.mp4", dpi=120):
    for f in files:
        hdu = fits.open(f)
        data = hdu[0].data
        header = hdu[0].header
        date_obs = header.get("DATE-OBS", "UNKNOWN")

        V = data[wave_V, 3, roi[0]:roi[1], roi[2]:roi[3]]*100

        dummy.set_data(V)
        ax.set_title(f"Stokes V [%] {date_obs}", fontsize=10)

        writer.grab_frame()
        hdu.close()

plt.close(fig)
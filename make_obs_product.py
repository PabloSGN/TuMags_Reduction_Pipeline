#!/usr/bin/env python3
# ./make_obs_products.py --fits "*Mg*.fits"
# ~/Nextcloud/PROGRAMACION_WORKINGFOLDER/TUMAG_DEVELOP/TuMags_Reduction_Pipeline/TuMags_Reduction_Pipeline/make_obs_product.py --fits "*Fe2.02*LV_1.0_v0.4*fits"

import os
import glob
import argparse
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import numpy as np
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams['animation.ffmpeg_path'] = '/opt/miniconda3/envs/torchmfbd-env/bin/ffmpeg'


# =========================
# ARGUMENTOS CLI
# =========================
parser = argparse.ArgumentParser(description="Generate products from FITS")

parser.add_argument(
    "--fits",
    required=True,
    help='Pattern for FITS files, e.g. "*Mg*.fits"'
)

parser.add_argument(
    "--out",
    default="output",
    help="Output directory"
)

args = parser.parse_args()

# =========================
# FILES
# =========================
files = sorted(glob.glob(args.fits))
files = [f for f in files if not os.path.basename(f).startswith("._")]

if len(files) == 0:
    raise RuntimeError("No FITS files found")

print(f"Found {len(files)} files")

# =========================
# AUTO-DETECCIÓN
# =========================
fname = os.path.basename(files[0])

# OBS_ID → antes de _TM_
obs_match = re.split(r"_TM_", fname)
obs_id = obs_match[0] if len(obs_match) > 1 else "UNKNOWN"

# VERSION → v...
ver_match = re.search(r"(v\d+\.?\d*)", fname)
version = ver_match.group(1) if ver_match else "v?"

# LINE → Mg / Fe...
if "Mg" in fname:
    line = "Mg"
elif "Fe" in fname:
    line = "Fe"
else:
    line = "UNKNOWN"

print("OBS_ID:", obs_id)
print("LINE:", line)
print("VERSION:", version)

# =========================
# CONFIG
# =========================
clim_V = (-1, 1)
clim_I = (0.6, 1.25)

roi = [0, -1, 0, -1]
wave_I = -1
wave_V = 2

# OUT DIRS
ref_dir = os.path.join(args.out, "reference_pngs")
mov_dir = os.path.join(args.out, "movies")

os.makedirs(ref_dir, exist_ok=True)
os.makedirs(mov_dir, exist_ok=True)

# =========================
# FUNCIONES
# =========================
def plot_clean(image, clim, title, cbar_label, outname,
               show_cbar=True, show_title=True):

    fig, ax = plt.subplots(figsize=(4, 4))

    im = ax.imshow(image, cmap="gray", clim=clim)
    ax.axis("off")

    if show_title:
        ax.set_title(title, fontsize=10)

    if show_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.02)
        cbar = plt.colorbar(im, cax=cax)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label(cbar_label, fontsize=9)

    plt.savefig(outname, dpi=150, bbox_inches="tight", pad_inches=0)
    plt.close()


def get_I_V(data):
    I = data[wave_I, 0, roi[0]:roi[1], roi[2]:roi[3]]
    V = data[wave_V, 3, roi[0]:roi[1], roi[2]:roi[3]] * 100
    return I, V


# =========================
# REFERENCIAS
# =========================
print("Generating reference images...")

f_first = files[0]
f_mid = files[len(files)//2]

# CONTINUUM
data = fits.getdata(f_first)
I, _ = get_I_V(data)

plot_clean(
    I, clim_I,
    f"{line} – Stokes I",
    "[Intensity]",
    f"{ref_dir}/{obs_id}_continuum.png"
)

# THUMBNAIL
plot_clean(
    I, clim_I,
    "",
    "",
    f"{ref_dir}/{obs_id}_continuum_thumb.png",
    show_cbar=False,
    show_title=False
)

# HIGHLIGHT
data = fits.getdata(f_mid)
_, V = get_I_V(data)

plot_clean(
    V, clim_V,
    f"{line} – Stokes V",
    "[%]",
    f"{ref_dir}/{obs_id}_highlight.png"
)

print("Generating details PNG...")

data = fits.getdata(files[0])
header = fits.getheader(files[0])
date_obs = header.get("DATE-OBS", "UNKNOWN")

I, V = get_I_V(data)

fig, axes = plt.subplots(1, 2, figsize=(8, 3))  # más compacto

for ax, img, clim, title in zip(
    axes,
    [I, V],
    [clim_I, clim_V],
    [f"{line} I", f"{line} V"]
):
    im = ax.imshow(img, cmap="gray", clim=clim)
    ax.set_title(title, fontsize=9)
    ax.axis("off")

plt.subplots_adjust(
    left=0.02,
    right=0.98,
    top=0.85,
    bottom=0.05,
    wspace=0.05   # 👈 separa muy poco los paneles
)

# título arriba sin estropear el layout
fig.text(0.5, 0.95, date_obs, ha='center', va='center', fontsize=10)

png_name = f"{obs_id}_TUMAG_details.png"

plt.savefig(
    os.path.join(mov_dir, png_name),
    dpi=150
)

plt.close()


# =========================
# MOVIE
# =========================
print("Generating movie...")

fps = 2

FFMpegWriter = manimation.writers['ffmpeg']
writer = FFMpegWriter(fps=fps)


# writer = FFMpegWriter(
#     fps=2,
#     codec='mpeg4',   # <- ESTE es el importante
#     bitrate=2000
# )

movie_name = f"{obs_id}_TUMAG_details.mp4"

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
with writer.saving(fig,os.path.join(mov_dir, movie_name), dpi=120):
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
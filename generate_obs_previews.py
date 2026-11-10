import os
import sys
import shutil
import argparse
import yaml
import numpy as np
from pathlib import Path
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ---------------------------------------------------------
# Configuración dinámica de ffmpeg
# ---------------------------------------------------------
# ---------------------------------------------------------
# Configuración estricta de ffmpeg (Prioriza Conda / sistema)
# ---------------------------------------------------------
conda_ffmpeg = os.path.join(sys.prefix, 'bin', 'ffmpeg')

if os.path.exists(conda_ffmpeg) and os.access(conda_ffmpeg, os.X_OK):
    # Fuerza el uso del binario del entorno activo (/opt/miniconda3/envs/torchmfbd-env/bin/ffmpeg)
    plt.rcParams['animation.ffmpeg_path'] = conda_ffmpeg
elif shutil.which('ffmpeg'):
    plt.rcParams['animation.ffmpeg_path'] = shutil.which('ffmpeg')
else:
    raise RuntimeError("No se encontró ejecutable de ffmpeg en el entorno de conda ni en el PATH.")

try:
    FFMpegWriter = manimation.writers['ffmpeg']
except KeyError:
    if shutil.which('ffmpeg'):
        plt.rcParams['animation.ffmpeg_path'] = shutil.which('ffmpeg')
        FFMpegWriter = manimation.writers['ffmpeg']
    else:
        raise RuntimeError("No se encontró 'ffmpeg'. Instálalo en el entorno con 'conda install -c conda-forge ffmpeg'.")


def load_yaml_config(yaml_file):
    with open(yaml_file, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    return config


def get_fits_files(data_dir, obs_id, line, proc_version, level_tag="LV_1.0"):
    """
    Busca y ordena los archivos FITS filtrando estrictamente por el nivel (level_tag).
    """
    if not os.path.exists(data_dir):
        data_dir = "./"

    all_files = sorted(os.listdir(data_dir))
    
    # Filtrado estricto por versión y etiqueta de nivel (ej: LV_1.0 o LV_1.2)
    filtered = [
        os.path.join(data_dir, f) for f in all_files
        if f.endswith('.fits') 
        and not f.startswith('._')
        and f'v{proc_version}' in f
        and line in f
        and level_tag in f
    ]
    
    # Si no encuentra con el tag exacto 'LV_1.0', busca alternativas tipo 'LV1.0' o '1.0'
    if not filtered:
        tag_clean = level_tag.replace("LV_", "").replace("LV", "")
        filtered = [
            os.path.join(data_dir, f) for f in all_files
            if f.endswith('.fits') 
            and not f.startswith('._')
            and f'v{proc_version}' in f
            and line in f
            and tag_clean in f
        ]
        
    if line == 'Mg' and len(filtered) > 1:
        filtered = filtered[1:]
        
    return sorted(filtered)


def extract_2d_slice(data, wave_idx, pol_idx, roi):
    """
    Extrae una sub-matriz bidimensional 2D limpia (Y, X) sin importar 
    cuántas dimensiones tenga el cubo FITS (4D o 5D).
    """
    # Si la estructura es 5D: (cam, wave, mod, y, x)
    if data.ndim == 5:
        slice_2d = data[0, wave_idx, pol_idx, roi[0]:roi[1], roi[2]:roi[3]]
    # Si la estructura es 4D: (wave, mod, y, x)
    elif data.ndim == 4:
        slice_2d = data[wave_idx, pol_idx, roi[0]:roi[1], roi[2]:roi[3]]
    # Si la estructura es 3D: (mod, y, x)
    elif data.ndim == 3:
        slice_2d = data[pol_idx, roi[0]:roi[1], roi[2]:roi[3]]
    else:
        slice_2d = data[roi[0]:roi[1], roi[2]:roi[3]]

    return np.squeeze(slice_2d)


def plot_single_stokes(image, clim, title, cbar_label, outname, date_obs, is_thumbnail=False):
    fig, ax = plt.subplots(figsize=(4, 4))
    im = ax.imshow(image, cmap="gray", clim=clim)
    ax.axis("off")

    if not is_thumbnail:
        ax.set_title(f"{title}\n{date_obs}", fontsize=10)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label(cbar_label, fontsize=9)

    plt.tight_layout(pad=0.1 if is_thumbnail else 1.08)
    plt.savefig(outname, dpi=150 if not is_thumbnail else 100, bbox_inches='tight' if is_thumbnail else None)
    plt.close(fig)


def generate_reference_pngs(files, obs_id, line, roi, clim_I, clim_V, wave_I, wave_V, out_dir="reference_pngs"):
    os.makedirs(out_dir, exist_ok=True)
    
    with fits.open(files[0]) as hdu:
        data = hdu[0].data
        header = hdu[0].header
        date_obs = header.get("DATE-OBS", "UNKNOWN")

        # wave_I = 0
        # wave_V = 2 if data.shape[0] > 2 else 0
        
        I_first = extract_2d_slice(data, wave_I, 0, roi)
        V_first = extract_2d_slice(data, wave_V, 3, roi) * 100.0

        # Auto-clim si es None
        c_I = clim_I if clim_I is not None else get_auto_clim(I_first, (0.5, 99.5))
        if clim_V is not None:
            c_V = clim_V
        else:
            v_max = np.nanpercentile(np.abs(V_first), 99)
            c_V = (-v_max, v_max)

        # Si Stokes V está centrado en cero, es mejor forzar una escala simétrica para el mapa de polarización:
        v_abs_max = np.nanpercentile(np.abs(V_first), 99)
        clim_V = (-v_abs_max, v_abs_max)

    plot_single_stokes(
        I_first, c_I,
        f"Stokes I ({line})", "[Intensity]",
        os.path.join(out_dir, f"{obs_id}_continuum.png"),
        date_obs
    )
    
    plot_single_stokes(
        I_first, c_I,
        "", "",
        os.path.join(out_dir, f"{obs_id}_thumbnail.png"),
        date_obs,
        is_thumbnail=True
    )

    with fits.open(files[len(files) // 2]) as hdu:
        data = hdu[0].data
        header = hdu[0].header
        date_obs = header.get("DATE-OBS", "UNKNOWN")
        
        V_mid = extract_2d_slice(data, wave_V, 3, roi) * 100.0

    plot_single_stokes( 
        V_mid, c_V,
        f"Stokes V Highlight ({line})", "[%]",
        os.path.join(out_dir, f"{obs_id}_highlight.png"),
        date_obs
    )

    print(f"[INFO] PNGs de referencia y Thumbnail guardados en: {out_dir}/")


def generate_stokes_movies(files, obs_id, line, roi, clim_I, clim_V, wave_I, wave_V, fps=2, out_dir="movies"):
    
    os.makedirs(out_dir, exist_ok=True)
    writer = FFMpegWriter(fps=fps)
    
    with fits.open(files[0]) as hdu0:
        data0 = hdu0[0].data
        init_I = extract_2d_slice(data0, wave_I, 0, roi)
        init_V = extract_2d_slice(data0, wave_V, 3, roi) * 100.0

    c_I = clim_I if clim_I is not None else get_auto_clim(init_I, (0.5, 99.5))
    
    if clim_V is not None:
        c_V = clim_V
    else:
        v_max = np.nanpercentile(np.abs(init_V), 99)
        c_V = (-v_max, v_max)

    # ---------- Película Stokes I ----------
    fig, ax = plt.subplots(figsize=(4, 4))
    dummy_I = ax.imshow(init_I, cmap="gray", clim=c_I)
    ax.axis("off")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(dummy_I, cax=cax)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label("[Ic]", fontsize=9)

    mov_path_I = os.path.join(out_dir, f"{obs_id}_{line}_details_StokesI.mp4")
    with writer.saving(fig, mov_path_I, dpi=120):
        for f in files:
            with fits.open(f) as hdu:
                data = hdu[0].data
                header = hdu[0].header
                date_obs = header.get("DATE-OBS", "UNKNOWN")
                I = extract_2d_slice(data, wave_I, 0, roi)
                dummy_I.set_data(I)
                ax.set_title(f"Stokes I ({line})\n{date_obs}", fontsize=10)
                writer.grab_frame()
    plt.close(fig)

    # ---------- Película Stokes V ----------
    fig, ax = plt.subplots(figsize=(4, 4))
    dummy_V = ax.imshow(init_V, cmap="gray", clim=c_V)
    ax.axis("off")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(dummy_V, cax=cax)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label("[%]", fontsize=9)

    mov_path_V = os.path.join(out_dir, f"{obs_id}_{line}_details_StokesV.mp4")
    with writer.saving(fig, mov_path_V, dpi=120):
        for f in files:
            with fits.open(f) as hdu:
                data = hdu[0].data
                header = hdu[0].header
                date_obs = header.get("DATE-OBS", "UNKNOWN")
                V = extract_2d_slice(data, wave_V, 3, roi) * 100.0
                dummy_V.set_data(V)
                ax.set_title(f"Stokes V [%] ({line})\n{date_obs}", fontsize=10)
                writer.grab_frame()
    plt.close(fig)

    print(f"[INFO] Animaciones MP4 guardadas en: {out_dir}/")

def get_auto_clim(image_data, percentiles=(1, 99)):
    """
    Calcula los límites vmin y vmax usando percentiles
    para evitar que los valores extremos oscurezcan el resto.
    """
    valid_data = image_data[~np.isnan(image_data)]
    if valid_data.size == 0:
        return (0, 1)
    vmin, vmax = np.nanpercentile(valid_data, percentiles)
    return (vmin, vmax)

def main():
    parser = argparse.ArgumentParser(description="Generador de visualizaciones para observaciones de TuMag.")
    
    parser.add_argument("yaml", nargs="?", default="config_tumag_13_PLAG_CLV_Mg_v0.46.yaml", help="Ruta al archivo YAML de configuración.")
    parser.add_argument("level", nargs="?", default="LV_1.0", help="Nivel FITS (ej. LV_1.0, LV_1.2).")
    
    parser.add_argument("--level", "-l", dest="level_opt", default=None, help="Nivel FITS alternativo en flag.")
    parser.add_argument("--wave-i", type=int, default=None, help="Índice de longitud de onda para Stokes I.")
    parser.add_argument("--wave-v", type=int, default=None, help="Índice de longitud de onda para Stokes V.")
    parser.add_argument("--clim-i", type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help="Límites de color para Stokes I: --clim-i 0.6 1.25")
    parser.add_argument("--clim-v", type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help="Límites de color para Stokes V: --clim-v -1.0 1.0")

    args = parser.parse_args()

    yaml_path = args.yaml
    level_tag = args.level_opt if args.level_opt is not None else args.level

    cfg = load_yaml_config(yaml_path)
    
    obs_id = cfg.get("obs_ID", "OBS_UNKNOWN")
    proc_version = cfg.get("proc_version", "0.46")
    line = str(cfg.get("process_line", "1"))
    
    if line == "1":
        line_label = "Mg"
    elif line == "2.02":
        line_label = "Fe2.02"
    else:
        line_label = line

    output_folder = cfg.get("output_folder", "./reduccion/")
    data_dir = os.path.join(output_folder, obs_id)
    roi = cfg.get("plots", {}).get("roi_plots", [0, -1, 0, -1])

    files = get_fits_files(data_dir, obs_id, line_label, proc_version, level_tag=level_tag)
    if not files:
        print(f"[ERROR] No se encontraron archivos FITS para {obs_id} en '{data_dir}' con v{proc_version}, línea '{line_label}' y nivel '{level_tag}'.")
        return

    # Determinar wave_I y wave_V por defecto según la forma del archivo si no se especificaron
    with fits.open(files[0]) as hdu0:
        shape = hdu0[0].data.shape
        wave_I = args.wave_i if args.wave_i is not None else (0 if line != 'Mg' else -1)
        wave_V = args.wave_v if args.wave_v is not None else (2 if shape[0] > 2 else 0)

    print(f"==================================================")
    print(f" PROCESANDO OBSERVACIÓN: {obs_id} (Línea: {line_label} | Nivel: {level_tag})")
    print(f" Parámetros: wave_I={wave_I}, wave_V={wave_V}, clim_I={args.clim_i}, clim_V={args.clim_v}")
    print(f"==================================================")
    print(f"[INFO] Se han encontrado {len(files)} archivos FITS para el nivel {level_tag}.")

    clim_I = tuple(args.clim_i) if args.clim_i else None
    clim_V = tuple(args.clim_v) if args.clim_v else None

    generate_reference_pngs(files, obs_id, line_label, roi, clim_I, clim_V, wave_I, wave_V)
    generate_stokes_movies(files, obs_id, line_label, roi, clim_I, clim_V, wave_I, wave_V)

if __name__ == "__main__":
    main()
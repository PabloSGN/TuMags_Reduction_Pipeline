#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TuMags Reduction Pipeline — main entry point (refactor)
"""

from __future__ import annotations

# ============================= STANDARD LIBS ================================= #
import os
import sys
import logging
import multiprocessing
from argparse import ArgumentParser, Namespace
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Optional

# ============================= THIRD-PARTY =================================== #
import numpy as np
import pandas as pd
from tqdm import tqdm
from astropy.io import fits

# ============================= PROJECT PATH ================================== #
def add_project_root() -> Path:
    """
    Inserta la carpeta del script en sys.path (raíz del proyecto).
    Evita rutas absolutas hardcodeadas locales.
    """
    this = Path(__file__).resolve()
    root = this.parent
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    return root

PROJECT_ROOT = add_project_root()

# =============================== LOCAL MODS ================================== #
# Carga de módulos locales (después de fijar PROJECT_ROOT)
import config as cf
import image_handler as ih
from master_dark import compute_master_darks
from master_flatfield import compute_master_flat_field
from image_filtering import filter_frecuencies
from fits_files_handling import generate_fits, update_header
from alignment import align_obsmode
from advanced_alignment import apply_transform, update_alignment_csv, run_alignment_for_wavelength
from demodulation import demodulate
from destretch import destretch
import pd_functions_v22 as phased

from xtalk_new import (
    fit_mueller_matrix, fit_mueller_matrix_tiled,
    write_crosstalk_header, fit_interference_Iref_to_Q_tiled,
    write_interference_Iref_to_Q_header, apply_crosstalk_coeffs_standard
)

from process_data_utils import (
    ConfigLoader, DEFAULT_CONFIG, parse_range,
    print_shifts_by_cam, plt_darks, plt_flats, plt_level,
    _timestamp_from_filename, balance, parse_header_time, minutes_from_dt,
    interpolate_filter, correct_image_fft, extract_crosstalk_coeffs
)

from process_data_timelines import obs_dict  # timeline_names y obs_dict

# ============================ THREADS & LOGGING ============================== #
def configure_threads(n: Optional[int] = None) -> None:
    """
    Limita hilos de BLAS para evitar oversubscription (útil con multiprocessing).
    Si 'n' es None, intenta respetar OMP_NUM_THREADS; por defecto 1.
    """
    val = str(n or os.environ.get("OMP_NUM_THREADS") or 1)
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(k, val)

def setup_logging(workspace: Path, level: int = logging.INFO) -> Path:
    """
    Configura logging con consola (INFO+) y fichero (WARNING+).
    Devuelve la ruta al fichero de log.
    """
    log_path = workspace / "errors.log"
    fmt = "%(asctime)s [%(processName)s] %(levelname)s: %(message)s"

    # Consola (root handler)
    logging.basicConfig(
        level=level,
        format=fmt,
        handlers=[logging.StreamHandler(sys.stdout)]
    )

    # Fichero (solo WARNING+)
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setLevel(logging.WARNING)
    fh.setFormatter(logging.Formatter(fmt))
    logging.getLogger().addHandler(fh)

    return log_path

# ============================ ARGUMENT PARSING =============================== #
def parse_args(argv=None) -> Namespace:
    p = ArgumentParser(
        description="TuMags Reduction Pipeline",
        epilog="Ejemplo: python3 process_data_main.py -f ./config_tumag.yaml"
    )
    p.add_argument(
        "-f", "--config",
        default="config_tumag.yaml",
        help="Ruta al archivo YAML de configuración"
    )
    p.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Número de hilos BLAS/numexpr (por defecto 1 si no se define)"
    )
    p.add_argument(
        "--no-ask-update",
        action="store_true",
        help="No preguntar para actualizar el YAML con defaults faltantes"
    )
    return p.parse_args(argv)

# ============================== GLOBAL STATE ================================= #
try:
    workspacePath = Path(__file__).resolve().parent
except Exception:
    workspacePath = Path("./").resolve()

errLogFilename = setup_logging(workspacePath, level=logging.INFO)
cam_linearity = ([1539, 1540], [1.0, 1.0])
# cam_linearity = 0

# ======================= processing programs ======================= #
# =======================  reduce_image_0_5   ======================= #
def reduce_image_0_5(ocs, OCs, cfg, dc_real, ff_data, obs_ID, ff_paths, dc_paths, process_line_index):
    process_name = multiprocessing.current_process().name
    logging.basicConfig(
        level=logging.INFO,
        format=f'%(asctime)s [{process_name}] %(levelname)s: %(message)s',
        force=True
    )

    logging.info(f' processing ocs: {ocs} ........... ')
    om_value = OCs[ocs]['OM']
    if cfg['process_line'] == om_value:
        logging.info(f' processing ocs: {ocs} which corresponds to {om_value} with {len(OCs[ocs]["ims"])} total images')
        if len(OCs[ocs]['ims']) != obs_dict[obs_ID]['obs_size'][process_line_index]:
            logging.error(f"  >> Error. The ocs: {ocs} number of images {len(OCs[ocs]['ims'])} does not coincide with the timeline info: {obs_dict[obs_ID]['obs_size'][process_line_index]}")
            return

        obs_data = ih.nominal_observation(cfg['process_line'], OCs[ocs]["ims"], dc_real, modify_linearity=cam_linearity)
        data = obs_data.get_data()
        om_info = obs_data.get_info()

        date_str = om_info["Images_headers"]["wv_0"]["M0"]["Date"].strftime("%d%m%YT%H%M%S")
        dataid = obs_ID + "_TM_" + cf.om_config[cfg['process_line']]["name"] + '_' + str(cf.om_config[cfg['process_line']]["Nlambda"]) + "_"
        filename = dataid + date_str
        extended_filename = f"{filename}_LV_0.5_v{cfg['proc_version']}.fits"
        logging.info(f' Output filename: {extended_filename}')

        logging.info(f'  >> FF correction..........')
        data = np.where(np.isfinite(ff_data), data / ff_data, 0)

        logging.info(f'  >> data cropping..........')
        data = data[:, :, :, cfg['centro'][0] - cfg['corte']:cfg['centro'][0] + cfg['corte'],
                    cfg['centro'][1] - cfg['corte']:cfg['centro'][1] + cfg['corte']]

        if cfg['level_05']['filtering']:
            data = filter_frecuencies(data, band='fixed', verbose=True, pad=500, N=350, cam=1)

        generate_fits(
            data.astype(np.float32),
            cfg['output_folder'] + obs_ID + '/', extended_filename, '0.5', cfg['proc_version'],
            om_info=om_info, zkes=None, shifts=None, fitted_muller=None, datatype='SCIENCE',
            DARK_ID=dc_paths, FLAT_ID=ff_paths[process_line_index]
        )

        if cfg['plots']['plot_level0_5']:
            plt_level(
                data, cfg['plots']['roi_plots'], cfg['output_folder'] + obs_ID,
                f"{filename}_LV_0.5_v{cfg['proc_version']}", '0.5', 'flat_corrected'
            )

# =======================  reduce_image_0_7   ======================= #
def reduce_image_0_7(input_data_filename, cfg, df, line, from_label="LV_0.5", to_label="LV_0.7"):
    # Subproceso: logging local
    process_name = multiprocessing.current_process().name
    logging.basicConfig(
        level=logging.INFO,
        format=f'%(asctime)s [{process_name}] %(levelname)s: %(message)s',
        force=True
    )

    logging.info(f' processing file: {input_data_filename} ')
    filename = Path(input_data_filename.replace(from_label, to_label)).stem
    obs_ID = cfg['obs_ID']

    # Leer data
    with fits.open(input_data_filename) as hdul:
        data = hdul[0].data
        header = hdul[0].header

    cn, wn, pn, xs, ys = data.shape
    roi = cfg['level_07']['align_roi']
    if roi[1] == -1 or roi[-1] == -1:
        roi[1] = ys
        roi[3] = xs

    scale, gamma = balance(data[0, 0], data[1, 0], roi=roi, clip_percentiles=(40, 60))
    data[1] = data[1] * scale
    logging.info(f"Balance (4x): scale={scale:.6f}, gamma={gamma:.6f}")

    timestamp = _timestamp_from_filename(input_data_filename)
    logging.info(f'  timestamp: {timestamp} ')

    # Debug opcional de filtros
    if cfg['debug']['filter_test']:
        cm, wv, pl, _, _ = data.shape
        data_ = np.copy(data)
        if line['line'] == '525.02':
            logging.info(f'  FILTERING NEW FREQ IN 525.02......')
            for cam in range(cm):
                for wl in range(wv):
                    for pn in range(pl):
                        if wl == 0:
                            kx = (22, 22, 22, 23, 23, 23, 23, 24, 24, 24)
                            ky = (-3, -4, -5, -2, -3, -4, -5, -2, -3, -4)
                            data_[cam, wl, pn, :, :] = correct_image_fft(
                                data[cam, wl, pn, :, :], kx=kx, ky=ky, h=1, mode='intensity'
                            )
                        if wl >= 0 and wl <= 6:
                            kx = (22, 22, 22, 23, 23, 23, 23, 24, 24, 24)
                            ky = (3, 4, 5, 2, 3, 4, 5, 2, 3, 4)
                            data_[cam, wl, pn, :, :] = correct_image_fft(
                                data[cam, wl, pn, :, :], kx=kx, ky=ky, h=1, mode='intensity'
                            )
                        if wl == 4:
                            kx = (16, 17, 18, 16, 17, 18)
                            ky = (-3, -3, -3, -4, -4, -4)
                            data_[cam, wl, pn, :, :] = correct_image_fft(
                                data[cam, wl, pn, :, :], kx=kx, ky=ky, h=1, mode='intensity'
                            )
                        if wl == 7:
                            kx = (20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24)
                            ky = (-4, -5, -6, -7, -4, -5, -6, -7, -3, -4, -5, -6, -7, -3, -4, -5, -6, -7, -3, -4, -5, -6, -7)
                            data_[cam, wl, pn, :, :] = correct_image_fft(
                                data[cam, wl, pn, :, :], kx=kx, ky=ky, h=1, mode='intensity'
                            )
        if line['line'] == '517':
            logging.info(f'  FILTERING NEW FREQ IN 517......')
            for cam in range(cm):
                for wl in range(wv):
                    for pn in range(pl):
                        if wl == 0:
                            kx = (22, 22, 22, 23, 23, 23, 24, 24, 24)
                            ky = (-2, -3, -4, -2, -3, -4, -2, -3, -4)
                            data_[cam, wl, pn, :, :] = correct_image_fft(
                                data[cam, wl, pn, :, :], kx=kx, ky=ky, h=1, mode='intensity'
                            )
                        if wl >= 7 and wl <= 8:
                            kx = (22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24)
                            ky = (-2, -3, -4, -2, -3, -4, -5, -2, -3, -4, -5)
                            data_[cam, wl, pn, :, :] = correct_image_fft(
                                data[cam, wl, pn, :, :], kx=kx, ky=ky, h=1, mode='intensity'
                            )
                        if wl == 9:
                            kx = (20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24)
                            ky = (-2, -3, -4, -5, -2, -3, -4, -5, -2, -3, -4, -5, -2, -3, -4, -5, -2, -3, -4, -5)
                            data_[cam, wl, pn, :, :] = correct_image_fft(
                                data[cam, wl, pn, :, :], kx=kx, ky=ky, h=1, mode='intensity'
                            )
        data = np.copy(data_)
        del data_

    def _print_result_pretty(result, title="Resultado interpolado"):
        print("\n" + "=" * 70)
        print(f" {title}")
        print("=" * 70)
        for key, value in result.items():
            print(f"{key:<12} : {value}")
        print("=" * 70 + "\n")

    # --- Normalizador DF (gestiona wave faltante como -1) ---
    def _normalize_alignment_df(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        required = [
            'timestamp',
            'angle_deg', 't_y', 't_x', 'center_y', 'center_x',
            'scale_x', 'scale_y', 'shear_x', 'shear_y',
            'fun', 'obs_ID', 'line', 'wave'
        ]
        for c in required:
            if c not in df.columns:
                df[c] = np.nan

        num_cols = [
            'timestamp', 'angle_deg', 't_y', 't_x', 'center_y', 'center_x',
            'scale_x', 'scale_y', 'shear_x', 'shear_y', 'fun'
        ]
        for c in num_cols:
            df[c] = pd.to_numeric(df[c], errors='coerce')

        df['obs_ID'] = df['obs_ID'].astype(str).fillna("")
        df['line'] = df['line'].astype(str).str.strip().fillna("")
        df['wave'] = pd.to_numeric(df.get('wave', -1), errors='coerce').fillna(-1).astype(int)

        df = df.sort_values(['obs_ID', 'line', 'wave', 'timestamp']).reset_index(drop=True)
        return df

    # =================== ADVANCED ALIGNMENT =================== #
    if cfg['level_07']['advanced_alignment']:
        logging.info("  Entering advanced_alignment .....")

        # Ruta del CSV
        MODULE_DIR = Path(__file__).resolve().parent
        alignment_csv_filename = MODULE_DIR / "CALDATA" / f"{cfg['level_07']['align_rot_data_filter']}"

        advanced_wave = cfg['level_07']['advanced_wave']
        line_value = line['line']
        wave_list = advanced_wave if isinstance(advanced_wave, (list, tuple)) else [advanced_wave]

        for wave in wave_list:
            logging.info(f"  Aligning advanced mode for λ={wave} line={line_value}")

            data, result, timestamp = run_alignment_for_wavelength(
                data, line_value, wave, header
            )

            logging.info(f"  Done advanced_alignment for λ={wave}")
            _print_result_pretty(result)

            if cfg['level_07']['advanced_save']:
                update_alignment_csv(
                    alignment_csv_filename,
                    timestamp,
                    result.x,
                    result.fun,
                    overwrite=cfg['level_07']['advanced_overwrite'],
                    obs_ID=cfg['obs_ID'],
                    line=line_value,
                    wave=wave,
                )

        # Recarga y normaliza CSV
        if alignment_csv_filename.exists():
            df = pd.read_csv(
                alignment_csv_filename,
                dtype={'obs_ID': 'string', 'line': 'string'},
                converters={'wave': lambda x: int(float(x)) if pd.notna(x) and str(x).strip() != '' else -1}
            )
            df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['timestamp'])
            df = _normalize_alignment_df(df)
            logging.info("  Reloaded & normalized alignment CSV after advanced mode.")
        else:
            logging.warning("  alignment_results.csv no existe aún tras advanced mode.")

    # =================== INTERPOLACIÓN (sin advanced) =================== #
    df = _normalize_alignment_df(df)

    # ---- Nivel: DATASET (una sola interpolación para todo el cubo) ---- #
    if cfg['level_07']['aligment_level'] == 'dataset':
        logging.info('  Interpolating rotation (dataset-level) ...')

        try:
            timestamp_ds = timestamp
        except NameError:
            dt = [parse_header_time(header[f"WV_{0}_M{k}"]) for k in range(4)]
            timestamp_ds = sum(minutes_from_dt(d) for d in dt) / 4.0

        result = interpolate_filter(
            df=df,
            timestamp=timestamp_ds,
            line=line['line'],
            obs_ID=cfg['obs_ID'],
            wave=None
        )
        _print_result_pretty(result)

        for i in range(wn):
            for j in range(pn):
                data[1, i, j] = apply_transform(
                    data[1, i, j],
                    result['angle_deg'],
                    np.array([result['t_y'], result['t_x']]),
                    np.array([result['center_y'], result['center_x']]),
                    scale_x=result["scale_x"], scale_y=result["scale_y"],
                    shear_x=result["shear_x"], shear_y=result["shear_y"]
                )

    # ---- Nivel: WAVE (interpolación por λ) ---- #
    elif cfg['level_07']['aligment_level'] == 'wave':
        logging.info("  Interpolating rotation (wave-level) ...")

        for i in range(wn):
            dt = [parse_header_time(header[f"WV_{i}_M{k}"]) for k in range(4)]
            timestamp_i = sum(minutes_from_dt(d) for d in dt) / 4.0

            result = interpolate_filter(
                df=df,
                timestamp=timestamp_i,
                line=line['line'],
                obs_ID=cfg['obs_ID'],
                wave=i  # si no hay, fallback a -1 dentro de interpolate_filter
            )
            _print_result_pretty(result)

            for j in range(pn):
                data[1, i, j] = apply_transform(
                    data[1, i, j],
                    result['angle_deg'],
                    np.array([result['t_y'], result['t_x']]),
                    np.array([result['center_y'], result['center_x']]),
                    scale_x=result["scale_x"], scale_y=result["scale_y"],
                    shear_x=result["shear_x"], shear_y=result["shear_y"]
                )

    # =================== Align mode (fourier/destretch) =================== #
    if cfg['level_07']['align_mode'] == 'fourier':
        try:
            verb = cfg['level_07']['align_verbose']
        except Exception:
            verb = False

        try:
            debug = cfg['debug']
        except Exception:
            debug = False

        try:
            align_modulations = cfg['level_07']['align_modulations']
        except Exception:
            align_modulations = True

        data, shifts, _ = align_obsmode(
            data, acc=cfg['level_07']['align_accuracy'], verbose=verb,
            filter=line['line'], returnshifts=True, roi=roi,
            quadrants=cfg['level_07']['align_quadrants'],
            align_sequence=cfg['level_07'].get('align_sequence', 0),
            debug=debug, align_modulations=align_modulations
        )

        if cfg['level_07']['align_quadrants'] == 0:
            update_header(header, 'REALIGN', 1)
            update_header(header, 'ALIGMETH', 'sicairos')
            update_header(header, 'ALIGN_ID', f"roi = {roi[0]}:{roi[1]},{roi[2]}{roi[3]}")

            for i in range(wn):
                for j in range(pn):
                    key = f'row0_{i}{j}'
                    val = float(shifts[i, 0, 0, j])
                    header[key] = (val, f'cam 0 row shift wn/pn [{i},{j}]')

                    key = f'col0_{i}{j}'
                    val = float(shifts[i, 0, 1, j])
                    header[key] = (val, f'cam 0 col shift wn/pn [{i},{j}]')

                    key = f'row1_{i}{j}'
                    val = float(shifts[i, 1, 0, j])
                    header[key] = (val, f'cam 1 row shift wn/pn [{i},{j}]')

                    key = f'col1_{i}{j}'
                    val = float(shifts[i, 1, 1, j])
                    header[key] = (val, f'cam 1 col shift wn/pn [{i},{j}]')

            print_shifts_by_cam(shifts)

    if cfg['level_07']['align_mode'] == 'destretch':
        data, _ = destretch(
            data,
            n_iterations=cfg['level_07']['destretch_n_iterations'],
            aling_cam=cfg['level_07']['destretch_aling_cam'],
            ngrid=cfg['level_07']['destretch_ngrid'],
            lr=cfg['level_07']['destretch_lr'],
            lambda_tt=cfg['level_07']['destretch_lambda_tt'],
            filter=line['line']
        )

    logging.info(f'  demodulation: ')
    # Modo dual I (opcional)
    if cfg['level_10']['crosst_dual']:
        _, data_both = demodulate(data, line['line'], dmod_matrices=cfg['level_07']['demod_matrix'], BothCams=True)
        out_file = input_data_filename.replace(from_label, to_label + '_stokesI' + cfg['level_07']['add_level_07_label'])
        logging.info(f' Saving filename: {out_file}')
        with fits.open(input_data_filename) as hdu_list:
            hdu_list[0].data = data_both[:, :, 0, :, :]
            hdu_list[0].header = header
            hdu_list.writeto(out_file, overwrite=True)

    # Demodulación principal
    data = demodulate(data, line['line'], dmod_matrices=cfg['level_07']['demod_matrix'], mode=cfg['level_07']['demod_mode'])

    if cfg['plots']['plot_level0_7']:
        plt_level(
            data, cfg['plots']['roi_plots'], cfg['output_folder'] + obs_ID,
            filename, '0.7', label='_' + cfg['level_07']['align_mode'] + 'demodulated' + cfg['level_07']['add_level_07_label']
        )

    out_file = input_data_filename.replace(from_label, to_label + cfg['level_07']['add_level_07_label'])
    logging.info(f' Saving filename: {out_file}')
    with fits.open(input_data_filename) as hdu_list:
        hdu_list[0].data = data
        hdu_list[0].header = header
        hdu_list.writeto(out_file, overwrite=True)

# =======================  reduce_image_1_0   ======================= #
def reduce_image_1_0(input_data_filename, cfg, from_label="LV_0.7", to_label="LV_1.0"):
    process_name = multiprocessing.current_process().name
    logging.basicConfig(
        level=logging.INFO,
        format=f'%(asctime)s [{process_name}] %(levelname)s: %(message)s',
        force=True
    )

    logging.info(f'  >> processing file: {input_data_filename} ')
    filename = Path(input_data_filename.replace(from_label, to_label)).stem
    obs_ID = cfg['obs_ID']

    with fits.open(input_data_filename) as hdul:
        data = hdul[0].data
        header = hdul[0].header

    if cfg['level_10']['crosst_dual']:
        filename_dual = Path(input_data_filename.replace('0.7', '0.7_stokesI')).with_suffix('')
        with fits.open(str(filename_dual) + '.fits') as hdul:
            dualI = hdul[0].data
    else:
        dualI = None

    roi = cfg['level_10']['crosst_roi']

    if cfg['level_10']['normalization'] == 0:
        norm_factor = np.median(data[-1, 0, roi[0]:roi[1], roi[2]:roi[3]])
        logging.info(f'  >> Normalization factor: {norm_factor} ')
    else:
        norm_factor = cfg['level_10']['normalization']
    data = data / norm_factor

    if cfg['level_10']['crosst_quadrants'] == 0:
        params, coeffs, per_wv = extract_crosstalk_coeffs(cfg['crosstalk'])

        if params['apply_only']:
            data, info_apply = apply_crosstalk_coeffs_standard(
                data, coeffs, per_wavelength=per_wv,
                use_local=params['use_local'], local_order=params['local_order'],
                deriv_sigma=params['deriv_sigma'], channels=params['channels']
            )
            logging.info(
                f"[crosstalk] apply_only={params['apply_only']}, use_local={params['use_local']}, "
                f"local_order={params['local_order']}, deriv_sigma={params['deriv_sigma']}, "
                f"channels={params['channels']}, per_wavelength={per_wv}"
            )
        else:
            pass

        data, info = fit_mueller_matrix(
            data,
            method=cfg['level_10']['crosst_mode'],
            pthresh=cfg['level_10']['crosst_threshold'],
            ctmethod='linfit',
            last_wvl=cfg['level_10']['crosst_last_wave'],
            region=cfg['level_10']['crosst_region'],
            use_local=cfg['level_10']['use_local'],
            local_order=1,
            deriv_sigma=0.0,
            local_ridge_lambda=0.0,
            aggregate_wavelengths=False,
            channels=('Q', 'U', 'V'),
            verbose=cfg['level_10']['crosst_verbose'],
            dualI=dualI,
        )

        if cfg['level_10']['crosst_mode'] == 'standard':
            header = write_crosstalk_header(header, info, channels=('Q', 'U', 'V'), index_width=3)

        if cfg['level_10']['crosst_mode'] == 'jaeggli':
            matrix = info.get('MM1a', {})
            update_header(header, 'CROSTALK', 1, after='ALIGMETH', comment='Was crosstalk correction applied?')
            for i in range(4):
                for j in range(4):
                    key = f'MMAT_{i}{j}'
                    val = float(matrix[i, j])
                    header[key] = (val, f'Mueller matrix element [{i},{j}]')

    else:
        data, info = fit_mueller_matrix_tiled(
            data,
            method=cfg['level_10']['crosst_mode'],
            pthresh=cfg['level_10']['crosst_threshold'],
            ctmethod='linfit',
            last_wvl=cfg['level_10']['crosst_last_wave'],
            region=cfg['level_10']['crosst_region'],
            use_local=cfg['level_10']['use_local'],
            local_order=1,
            deriv_sigma=0.0,
            local_ridge_lambda=0.0,
            aggregate_wavelengths=False,
            channels=('Q', 'U', 'V'),
            verbose=cfg['level_10']['crosst_verbose'],
            divisions=cfg['level_10']['crosst_quadrants'],
        )

    if cfg['level_10']['crosst_interference']:
        ref_wvl = -1
        data, tiles, out = fit_interference_Iref_to_Q_tiled(
            data,
            ref_wvl=ref_wvl,
            divisions=int(cfg['level_10']['crosst_interference']),
            region=[0, -1, 0, -1],
            ctmethod='linfit',
            n_sigma=3,
            pthresh_intensity=0.0,
            use_local=False,
            local_order=1,
            deriv_sigma=0.0,
            local_ridge_lambda=0.0,
            apply=True,
            show_grid=False,
            verbose=False
        )

        header = write_interference_Iref_to_Q_header(
            header,
            a=out['a'], b=out['b'], c=out['c'], d=out['d'], e=out['e'],
            ref_wvl=out['ref_wvl'], divisions=out['divisions'], tiles=tiles,
            use_local=out['use_local'], local_order=out['local_order'],
            deriv_sigma=out['deriv_sigma'], ridge_lambda=out['ridge'],
            index_width=3, after_keyword='ALIGMETH'
        )

    if cfg['plots']['plot_level1_0']:
        plt_level(
            data, cfg['plots']['roi_plots'], cfg['output_folder'] + obs_ID,
            filename, '1.0', 'demod' + cfg['level_10']['add_level_10_label']
        )

    out_file = input_data_filename.replace(from_label, to_label + cfg['level_10']['add_level_10_label'])
    logging.info(f' Saving filename: {out_file}')
    with fits.open(input_data_filename) as hdu_list:
        hdu_list[0].data = data
        hdu_list[0].header = header
        hdu_list.writeto(out_file, overwrite=True)

# =======================  reduce_image_1_1   ======================= #
def reduce_image_1_1(input_data_filename, cfg, zk=None):
    process_name = multiprocessing.current_process().name
    logging.basicConfig(
        level=logging.INFO,
        format=f'%(asctime)s [{process_name}] %(levelname)s: %(message)s',
        force=True
    )

    logging.info(f'  >> processing file: {input_data_filename} ')
    filename = Path(input_data_filename.replace("LV_1.0", "LV_1.1")).stem
    obs_ID = cfg['obs_ID']

    with fits.open(input_data_filename) as hdul:
        data = hdul[0].data
        header = hdul[0].header

    wn, pn, xs, ys = data.shape

    with tqdm(total=wn * pn) as pbar:
        for wl in range(wn):
            for pl in range(pn):
                if pl == 0:
                    data[wl, pl], noise_filter = phased.restore_ima(
                        data[wl, pl], zk, pd=0, low_f=0.2, reg1=0.05, reg2=1,
                        cobs=32.4, epsilon=0.02, sigma=5000, stray='moffat'
                    )
                else:
                    data[wl, pl], _ = phased.restore_ima(
                        data[wl, pl], zk, pd=0, low_f=0.2,
                        noise=noise_filter, reg1=0.05, reg2=1,
                        cobs=32.4, epsilon=0.02, sigma=5000, stray='moffat'
                    )
            pbar.update(1)

    if cfg['plots']['plot_level1_1']:
        plt_level(
            data, cfg['plots']['roi_plots'], cfg['output_folder'] + obs_ID,
            filename, '1.0', 'pd'
        )

    with fits.open(input_data_filename) as hdu_list:
        hdu_list[0].data = data
        hdu_list[0].header = header
        hdu_list.writeto(input_data_filename.replace("LV_1.0", "LV_1.1"), overwrite=True)

# =======================  reduce_image_0_6   ======================= #
def reduce_image_0_6(input_data_filename, cfg, zk=None, from_label="LV_0.5", to_label="LV_0.6"):
    process_name = multiprocessing.current_process().name
    logging.basicConfig(
        level=logging.INFO,
        format=f'%(asctime)s [{process_name}] %(levelname)s: %(message)s',
        force=True
    )

    logging.info(f'  >> processing file: {input_data_filename} ')
    filename = Path(input_data_filename.replace(from_label, to_label)).stem
    obs_ID = cfg['obs_ID']

    with fits.open(input_data_filename) as hdul:
        data = hdul[0].data
        header = hdul[0].header

    cam, wn, pn, xs, ys = data.shape

    zk_cam0 = np.copy(zk)
    zk_cam1 = np.copy(zk)
    zknew = [zk_cam0, zk_cam1]
    epsilon = [0.02, 0.02]

    with tqdm(total=wn * pn * cam) as pbar:
        for cm in range(cam):
            for wl in range(wn):
                for pl in range(pn):
                    if pl == 0:
                        data[cm, wl, pl], noise_filter = phased.restore_ima(
                            data[cm, wl, pl], zknew[cm], pd=0, low_f=0.2,
                            reg1=0.05, reg2=1, cobs=32.4, epsilon=epsilon[cm], sigma=5000, stray='moffat'
                        )
                    else:
                        data[cm, wl, pl], _ = phased.restore_ima(
                            data[cm, wl, pl], zknew[cm], pd=0, low_f=0.2, noise=noise_filter,
                            reg1=0.05, reg2=1, cobs=32.4, epsilon=epsilon[cm], sigma=5000, stray='moffat'
                        )
                    pbar.update(1)

    if cfg['plots']['plot_level0_6']:
        plt_level(
            data, cfg['plots']['roi_plots'], cfg['output_folder'] + obs_ID,
            f"{filename}_{from_label}_v{cfg['proc_version']}", '0.5', 'pd'
        )

    with fits.open(input_data_filename) as hdu_list:
        hdu_list[0].data = data
        hdu_list[0].header = header
        hdu_list.writeto(input_data_filename.replace(from_label, to_label), overwrite=True)

# ===================================================================
# =======================     MAIN PROGRAM    =======================
# ===================================================================
def main(argv=None) -> int:
    logging.info('-----------------------------------')
    logging.info('  >> Running process_data_main ')

    # --- CLI ---
    args = parse_args(argv)

    # --- Threads BLAS ---
    configure_threads(args.threads)

    # --- Multiprocessing (spawn) ---
    try:
        multiprocessing.set_start_method("spawn", force=True)
    except RuntimeError:
        pass

    # --- Cargar configuración con defaults y opción de actualizar YAML ---
    cfg_loader = ConfigLoader(
        filepath=args.config,
        defaults=DEFAULT_CONFIG,
    )
    cfg = cfg_loader
    cfg_dict = cfg.config

    obs_ID = cfg.obs_ID
    line = cfg.process_line
    logging.info(f"  >> Input obs_ID: {obs_ID}")
    logging.info(f"  >> Processing line obs_ID: {line}")
    logging.info(f"  >> Folder data location: {cfg.tumag_data_location}")
    logging.info(f"  >> Output folder location: {cfg.output_folder}")

    # ===================== LEVEL 0.5 / FLATS / DARKS ===================== #
    if cfg.force_redo['level0_5'] or cfg.force_redo['redo_flat'] or cfg.force_redo['redo_dark']:

        # Rutas base
        dc_paths = obs_dict[obs_ID]['darks'][0]
        ff_paths = obs_dict[obs_ID]['flats']
        obs_paths = obs_dict[obs_ID]['obsdata'][0]

        logging.info(f"  >> Darks: {dc_paths}, flats: {ff_paths}, obs: {obs_paths} ")
        logging.info('-----------------------------------')
        logging.info(f'  >> checking if {cfg.Organized_files_local_folder_name} exist and continuing')
        if not Path(cfg.tumag_data_location + cfg.Organized_files_local_folder_name).exists():
            logging.error(f' >> the Organized_files_local_folder_name: {cfg.Organized_files_local_folder_name} not found')
            return 1

        ih.Organization_folder_files = os.path.join(cfg.tumag_data_location, cfg.Organized_files_local_folder_name)

        logging.info(f'  >> checking if out folder {cfg.output_folder + obs_ID} exist and creating it if it does not')
        out_dir = Path(cfg.output_folder + obs_ID)
        if not out_dir.exists():
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / 'pngs').mkdir(parents=True, exist_ok=True)
        logging.info('-----------------------------------')

        # DARKS
        logging.info('  >> Running dark calculation ')
        dark_output_file = out_dir / f"{dc_paths}.npz"
        logging.info(f'  >> dark output file will be: {dark_output_file}')

        if dark_output_file.exists():
            logging.info(f'  >> dark output file exist.')
            if not cfg.force_redo["redo_dark"]:
                logging.info(f'  >> loading dark')
                dark = np.load(dark_output_file, allow_pickle=True)
                dc_real = dark['dc_real']
            else:
                logging.info(f'  >> but force_redo is True.')
                dark_paths = ih.get_images_paths(dc_paths)
                dc_real, head, rms_dark = compute_master_darks(dark_paths[cfg.darks_indexes["dark_from"]:cfg.darks_indexes["dark_to"]], verbose=True)
                np.savez(dark_output_file, dc_real=dc_real.astype(np.float32))
        else:
            logging.info(f' >> dark output file does not exist.')
            dark_paths = ih.get_images_paths(dc_paths)
            dc_real, head, rms_dark = compute_master_darks(dark_paths[cfg.darks_indexes["dark_from"]:cfg.darks_indexes["dark_to"]], verbose=True)
            np.savez(dark_output_file, dc_real=dc_real.astype(np.float32))

        if cfg.plots["plot_darks"]:
            plt_darks(dc_real)

        # FLATS
        logging.info('-----------------------------------')
        logging.info('  >> Running FLAT calculation ')

        # Índice de process_line
        try:
            obs_is_list = obs_dict[obs_ID]['obs_is']
            if cfg.process_line in obs_is_list:
                process_line_index = obs_is_list.index(cfg.process_line)
                logging.info(f"  >> process_line '{cfg.process_line}' found at position {process_line_index} in obs_is list.")
            else:
                logging.error(f"  >> process_line '{cfg.process_line}' not found in obs_is list: {obs_is_list}")
                return 1
        except Exception as e:
            logging.error(f"  >> Error checking process_line in obs_is list: {e}")
            return 1

        if cfg.flat_file:
            logging.info(f'  >> reading flat: {cfg.flat_file}')
            flat = np.load(cfg.flat_file, allow_pickle=True)
            ff_data = flat['ff_data']
            ff_info = flat['ff_info']
        else:
            flat_output_file = out_dir / f"{ff_paths[process_line_index]}.npz"
            logging.info(f'  >> flat output file will be: {flat_output_file}')

            if flat_output_file.exists():
                logging.info(f'  >> flat output file exist.')
                if not cfg.force_redo["redo_flat"]:
                    logging.info(f'  >> loading flat')
                    flat = np.load(flat_output_file, allow_pickle=True)
                    ff_data = flat['ff_data']
                    ff_info = flat['ff_info']
                else:
                    logging.info(f'  >> but force_redo is True.')
                    flat_paths = ih.get_images_paths(ff_paths[process_line_index])
                    ff_data, ff_info = compute_master_flat_field(
                        flat_paths, dc=dc_real, verbose=True,
                        modify_linearity=cam_linearity, norm_roi=cfg.flat_norm_roi,
                        norm_method=cfg.norm_method, remove_prefilter=cfg.remove_prefilter,
                        pref_model=cfg.pref_model, import_blueshift_guess=cfg.import_blueshift_guess
                    )
                    np.savez(flat_output_file, ff_data=ff_data.astype(np.float32), ff_info=ff_info)
            else:
                logging.info(f' >> flat output file does not exist.')
                flat_paths = ih.get_images_paths(ff_paths[process_line_index])
                ff_data, ff_info = compute_master_flat_field(
                    flat_paths, dc=dc_real, verbose=True,
                    modify_linearity=cam_linearity, norm_roi=cfg.flat_norm_roi,
                    norm_method=cfg.norm_method, remove_prefilter=cfg.remove_prefilter,
                    pref_model=cfg.pref_model, import_blueshift_guess=cfg.import_blueshift_guess
                )
                np.savez(flat_output_file, ff_data=ff_data.astype(np.float32), ff_info=ff_info)

            if cfg.plots["plot_flats"]:
                plt_flats(ff_data)

        logging.info('-----------------------------------')

        # OCs
        if cfg.force_redo['level0_5']:
            ocs_output_file = out_dir / f"{obs_ID}_ocs.npz"
            logging.info(f'  >> ocs output file will be: {ocs_output_file}')

            if ocs_output_file.exists():
                logging.info(f'  >> ocs output file exist.')
                OCs = np.load(ocs_output_file, allow_pickle=True)['OCs'][()]
            else:
                logging.info(f'  >> determining ocs....')
                obs_images = ih.get_images_paths(obs_paths)
                OCs = ih.separate_ocs(obs_images, verbose=False)
                np.savez(ocs_output_file, OCs=OCs)

            process_ocs = list(OCs.keys())
            logging.info(f'  >> available ocs {len(process_ocs)}')
            if len(process_ocs) == 0:
                logging.error(f'  >> Error. No ocs found for obs_ID: {obs_ID}')
                return 1

            process_ocs = [process_ocs[i] for i in parse_range(cfg.process_ocs, max_value=len(process_ocs)-1)]
            logging.info(f'  >> process ocs {process_ocs}')

            reduce_partial = partial(
                reduce_image_0_5,
                OCs=OCs, cfg=cfg_dict, dc_real=dc_real, ff_data=ff_data,
                obs_ID=obs_ID, ff_paths=ff_paths, dc_paths=dc_paths,
                process_line_index=process_line_index
            )

            if len(process_ocs) > 1 and cfg.parallel:
                with ProcessPoolExecutor(max_workers=cfg.max_workers) as executor:
                    executor.map(reduce_partial, process_ocs)
            elif len(process_ocs) > 1 and not cfg.parallel:
                for oc in process_ocs:
                    reduce_partial(oc)
            else:
                reduce_partial(process_ocs[0])

    logging.info('-----------------------------------')

    # LEVEL 0.6 (opcional PD antes)
    try:
        deconvolve_first = cfg.force_redo['level0_6']
        from_label = "LV_0.5"
        to_label = "LV_0.6"
    except Exception:
        deconvolve_first = False

    if deconvolve_first:
        logging.info(f'  >> Procesing level 0.6')
        cfg.force_redo['level1_1'] = False

        dataid = obs_ID + "_TM_" + cf.om_config[cfg.process_line]["name"] + '_' + str(cf.om_config[cfg.process_line]["Nlambda"]) + "_"
        ext_ = f"_{from_label}_v{cfg.proc_version}.fits"

        directory = Path(cfg.output_folder + obs_ID + '/')
        files_list = sorted([f for f in os.listdir(directory) if f.endswith('.fits') and not f.startswith('._')])
        files_list = [f for f in files_list if dataid in f and ext_ in f]
        files_list = [str(directory / f) for f in files_list]

        indices = parse_range(cfg.process_files, max_value=len(files_list)-1)
        if len(files_list) >= 1 and len(indices) <= len(files_list):
            files = [files_list[i] for i in indices]
            logging.info(f'  >> available files {len(files)}')
        if len(files_list) == 0:
            logging.error(f'  >> Error. No files found for obs_ID: {obs_ID}')
            return 1
        if len(indices) > len(files_list):
            logging.error(f'  >> More indices than files {indices} len of file list {len(files_list)}')
            return 1

        if not isinstance(files, list):
            files = [files]

        logging.info('Importing zernikes: %s', cfg.zernike_id)
        BASE_DIR = Path(__file__).resolve().parent
        file_path_csv = BASE_DIR / "TuMag_PD_results_All_filters_clean.csv"
        zk = phased.import_zernikes(cfg.zernike_id, csv_path=file_path_csv)

        reduce_partial = partial(
            reduce_image_0_6, cfg=cfg_dict, zk=zk,
            from_label=from_label, to_label=to_label
        )

        if len(files) > 1 and cfg.parallel:
            with ProcessPoolExecutor(max_workers=cfg.max_workers) as executor:
                executor.map(reduce_partial, files)
        elif len(files) > 1 and not cfg.parallel:
            for f in files:
                reduce_partial(f)
        else:
            reduce_partial(files[0])

    logging.info('-----------------------------------')

    # LEVEL 0.7 (alignment + demodulation)
    try:
        if cfg.force_redo['use_pd']:
            cfg.force_redo['level1_1'] = False
            from_label = "LV_0.6"
            to_label = "LV_0.8"
        else:
            from_label = "LV_0.5"
            to_label = "LV_0.7"
    except Exception:
        from_label = "LV_0.5"
        to_label = "LV_0.7"

    if cfg.force_redo['level0_7']:
        logging.info(f'  >> processing level 0.7 (alignment and demodulation)')

        dataid = obs_ID + "_TM_" + cf.om_config[cfg.process_line]["name"] + '_' + str(cf.om_config[cfg.process_line]["Nlambda"]) + "_"
        ext_ = f"_{from_label}_v{cfg.proc_version}.fits"

        directory = Path(cfg.output_folder + obs_ID + '/')
        files_list = sorted([f for f in os.listdir(directory) if f.endswith('.fits') and not f.startswith('._')])
        files_list = [f for f in files_list if dataid in f and ext_ in f]
        files_list = [str(directory / f) for f in files_list]

        indices = parse_range(cfg.process_files, max_value=len(files_list)-1)
        if len(files_list) >= 1 and len(indices) <= len(files_list):
            files = [files_list[i] for i in indices]
            logging.info(f'  >> available files {len(files)}')
        if len(files_list) == 0:
            logging.error(f'  >> Error. No files found for obs_ID: {obs_ID}')
            return 1
        if len(indices) > len(files_list):
            logging.error(f'  >> More indices than files {indices} len of file list {len(files_list)}')
            return 1

        BASE_DIR = Path(__file__).resolve().parent
        file_path_csv = BASE_DIR / "CALDATA" / f"{cfg.level_07['align_rot_data_filter']}"

        # Lectura robusta del CSV: line como string, wave con fallback -1
        df = pd.read_csv(
            file_path_csv,
            dtype={'obs_ID': 'string', 'line': 'string'},
            converters={'wave': lambda x: int(float(x)) if pd.notna(x) and str(x).strip() != '' else -1}
        )
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['timestamp'])
        df = df.sort_values('timestamp').reset_index(drop=True)

        if not isinstance(files, list):
            files = [files]

        reduce_partial = partial(
            reduce_image_0_7,
            cfg=cfg_dict, df=df, line=cf.om_config[cfg.process_line],
            from_label=from_label, to_label=to_label
        )

        if len(files) > 1 and cfg.parallel:
            with ProcessPoolExecutor(max_workers=cfg.max_workers) as executor:
                executor.map(reduce_partial, files)
        elif len(files) > 1 and not cfg.parallel:
            for f in files:
                reduce_partial(f)
        else:
            reduce_partial(files[0])

    logging.info('-----------------------------------')

    # LEVEL 1.0
    if cfg.force_redo['level1_0']:
        try:
            if cfg.force_redo['use_pd']:
                cfg.force_redo['level1_1'] = False
                from_label = "LV_0.8"
                to_label = "LV_1.2"
            else:
                from_label = "LV_0.7"
                to_label = "LV_1.0"
        except Exception:
            from_label = "LV_0.7"
            to_label = "LV_1.0"

        logging.info(f'  >> processing level 1.0')

        dataid = obs_ID + "_TM_" + cf.om_config[cfg.process_line]["name"] + '_' + str(cf.om_config[cfg.process_line]["Nlambda"]) + "_"
        ext_ = f"_{from_label}{cfg.level_07['add_level_07_label']}_v{cfg.proc_version}.fits"

        directory = Path(cfg.output_folder + obs_ID + '/')
        files = sorted([f for f in os.listdir(directory) if f.endswith('.fits') and not f.startswith('._')])
        files = [f for f in files if dataid in f and ext_ in f]
        files = [str(directory / f) for f in files]

        if len(files) > 1:
            files = [files[i] for i in parse_range(cfg.process_files, max_value=len(files)-1)]
        logging.info(f'  >> available files {len(files)}')
        if len(files) == 0:
            logging.error(f'  >> Error. No files found for obs_ID: {obs_ID}')
            return 1

        if not isinstance(files, list):
            files = [files]

        reduce_partial = partial(
            reduce_image_1_0, cfg=cfg_dict, from_label=from_label, to_label=to_label
        )

        if len(files) > 1 and cfg.parallel:
            with ProcessPoolExecutor(max_workers=cfg.max_workers) as executor:
                executor.map(reduce_partial, files)
        elif len(files) > 1 and not cfg.parallel:
            for f in files:
                reduce_partial(f)
        else:
            reduce_partial(files[0])

    logging.info('-----------------------------------')

    # LEVEL 1.1
    if cfg.force_redo['level1_1']:
        logging.info(f'  >> processing level 1.1')

        dataid = obs_ID + "_TM_" + cf.om_config[cfg.process_line]["name"] + '_' + str(cf.om_config[cfg.process_line]["Nlambda"]) + "_"
        ext_ = f"_LV_1.0{cfg.level_07['add_level_07_label']}_v{cfg.proc_version}.fits"

        directory = Path(cfg.output_folder + obs_ID + '/')
        files = sorted([f for f in os.listdir(directory) if f.endswith('.fits') and not f.startswith('._')])
        files = [f for f in files if dataid in f and ext_ in f]
        files = [str(directory / f) for f in files]

        if len(files) > 1:
            files = [files[i] for i in parse_range(cfg.process_files, max_value=len(files)-1)]
        logging.info(f'  >> available files {len(files)}')
        if len(files) == 0:
            logging.error(f'  >> Error. No files found for obs_ID: {obs_ID}')
            return 1

        if not isinstance(files, list):
            files = [files]

        logging.info('Importing zernikes: %s', cfg.zernike_id)
        BASE_DIR = Path(__file__).resolve().parent
        file_path_csv = BASE_DIR / "TuMag_PD_results_All_filters_clean.csv"
        zk = phased.import_zernikes(cfg.zernike_id, csv_path=file_path_csv)

        reduce_partial = partial(
            reduce_image_1_1, cfg=cfg_dict, zk=zk
        )

        if len(files) > 1 and cfg.parallel:
            with ProcessPoolExecutor(max_workers=cfg.max_workers) as executor:
                executor.map(reduce_partial, files)
        elif len(files) > 1 and not cfg.parallel:
            for f in files:
                reduce_partial(f)
        else:
            reduce_partial(files[0])

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        logging.exception("Fallo no controlado en process_data_main")
        sys.exit(1)
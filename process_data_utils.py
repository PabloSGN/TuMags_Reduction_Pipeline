import yaml
import logging
import time
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime 
import re 
import numpy as np
import logging
from scipy.ndimage import gaussian_filter
import copy
from pathlib import Path
import pandas as pd
import logging
from decimal import Decimal, InvalidOperation
from typing import Any


DEFAULT_CONFIG = {
    # ----------------------------------------------------------
    # GLOBAL PARAMETERS
    # ----------------------------------------------------------
    "obs_ID": "06_SPOT_AR1",
    "proc_version": "0.1",  # con cross-talk standard

    # Observation mode to process (en tu YAML lo dejas como string)
    "process_line": "1",

    # Rango de OCS / files: ":" -> todos (o rango "a:b" o índice único)
    "process_ocs": ":",
    "process_files": ":",

    "parallel": False,
    "max_workers": 4,

    # Rutas
    "tumag_data_location": "./",
    "Organized_files_local_folder_name": "Organized_files_local",
    "output_folder": "./reduccion/",

    # ----------------------------------------------------------
    # REPROCESSING OPTIONS
    # ----------------------------------------------------------
    "force_redo": {
        "redo_flat": False,
        "redo_dark": False,
        "level0_5": False,
        "level0_6": False,
        "level0_7": True,
        "level1_0": True,
        "level1_1": False,   # if phase_diversity = 'before' do nothing
        "use_pd": False      # needs level0_6: true before (or run before)
    },

    # ----------------------------------------------------------
    # PLOTTING OPTIONS
    # ----------------------------------------------------------
    "plots": {
        "plot_darks": False,
        "plot_flats": False,
        "plot_level0_5": True,
        "plot_level0_6": True,
        "plot_level0_7": False,
        "plot_level1_0": True,
        "plot_level1_1": True,
        "roi_plots": [200, -200, 200, -200],
    },

    # ----------------------------------------------------------
    # DARK AND FLAT OPTIONS
    # ----------------------------------------------------------
    "darks_indexes": {
        "dark_from": None,
        "dark_to": None,
    },

    # sugerido por Fbailen (toda la imagen, pero sin mascara)
    "flat_norm_roi": [250, -250, 250, -250],

    # archivos opcionales (NPY) si quieres forzar dark/flat
    "dark_file": False,
    "flat_file": False,

    "import_blueshift_guess": True,
    "norm_method": "blueshift",
    "remove_prefilter": True,
    "pref_model": True,

    # ----------------------------------------------------------
    # LEVEL 0.5 OPTIONS
    # ----------------------------------------------------------
    "level_05": {
        "filtering": True,       # Remove high-frequency content
        "save_level_05": True,   # Save processed level 0.5 data
    },

    # Parámetros generales de geometría
    "size": 2016,
    "centro": [973, 996],       # [size // 2 - 35, size // 2 - 12]
    "corte": 800,

    # ----------------------------------------------------------
    # LEVEL 0.7 OPTIONS
    # ----------------------------------------------------------
    "level_07": {
        "align_mode": "fourier",           # "fourier" or "destretch"
        "align_sequence": 1,               # 0 o 1 según comentario
        "align_roi": [0, -1, 0, -1],
        "align_accuracy": 0.001,
        "align_rot_data_filter": "alignment_results.csv",  # e.g. "pinholes_Fe2.02.csv"
        "align_quadrants": 0,              # -8 / -16 según versiones
        "align_verbose": False,
        "align_modulations": False,

        # Advanced alignment
        "advanced_alignment": True,
        "advanced_wave": 0,
        "advanced_save": True,
        "advanced_overwrite": True,

        # Nivel de aplicación del alineamiento
        "aligment_level": "dataset",          # "dataset", "wave" o "pol"

        # Destretch options
        "destretch_aling_cam": "same",     # 'all' | 'full' | 'same'
        "destretch_n_iterations": 100,
        "destretch_ngrid": 4,
        "destretch_lr": 0.1,
        "destretch_lambda_tt": 0.02,

        # Demodulación
        "demod_matrix": "demod_matrices_david",  # 'demod_matrices_acampos', 'demod_matrices', 'demod_matrices_david_ct'
        "demod_mode": "join",

        # Etiqueta añadida a resultados de level_07
        "add_level_07_label": "",
    },

    # ----------------------------------------------------------
    # LEVEL 1.0 OPTIONS
    # ----------------------------------------------------------
    "level_10": {
        "crosst_mode": "standard",              # "jaeggli" or "standard"
        "crosst_strategy": "sequential",
        "crosst_roi": [0, -1, 0, -1],    # Normalization region
        "crosst_region": [0, -1, 0, -1], # Crosstalk calculation region
        "crosst_threshold": 0.01,               # puede ser 1 float o lista para Q,U,V
        "crosst_intensity_threshold": 0.1,      # 0 evita usar umbra; 0 => no usar
        "crosst_last_wave": None,               # usar -1 para standard y -2 para jaeggli (según nota)
        # "use_local": True,
        "crosst_quadrants": 0,
        "crosst_dual": False,
        "crosst_interference": 0,               # o fichero 'comm1_cal.npz'
        "normalization": 0,
        "plot_crosst_method": False,
        "crosst_verbose": False,
        "add_level_10_label": "",               # cadena vacía por defecto
    },

    # Bloque de crosstalk explícito (aplicar/usar coeficientes)
    "crosstalk": {
        "apply_only": False,     # activa modo "aplicar sin estimar"
        "use_local": True,       # usa derivadas (Ix,Iy[,L]) al aplicar
        "local_order": 1,        # 1 -> Ix,Iy; 2 -> + Laplaciano
        "deriv_sigma": 0.0,      # suavizado de I para derivadas
        "coeffs_global": {
            "Q": {"a": 10.5, "b": 2.1e-3, "c": -1.2e-4, "d": 9.5e-5, "e": 0.0},
            "U": {"a":  0.0, "b": 1.8e-3, "c":  0.0,    "d": 0.0,    "e": 0.0},
            "V": {"a":  0.0, "b": 0.0,    "c":  0.0,    "d": 0.0,    "e": 0.0},
        },
        "channels": ["Q", "U", "V"],
        "aggregate_wavelengths": False,
        },

    # ----------------------------------------------------------
    # LEVEL 1.1 OPTIONS
    # ----------------------------------------------------------
    "zernike_id": "06_SPOT_Mg1_0",  # "06_SPOT_Fe2.02_0"

    # ----- DEBUGGING OPTIONS
    "debug": {
        "activate": False,
        "zoom_center": [600, 600],
        "zoom_size": 600,
        "filter_test": False,
    },
    "pd": 0,
    "low_f": 0.2,
    "reg1": 0.05,
    "reg2": 1,
    "cobs": 32.4,
    "epsilon": 0.02,
    "sigma": 5000,
    "stray": 'moffat'
}

def _pretty_value(v: Any) -> str:
    """
    Devuelve una representación legible del valor por defecto:
    - Para escalares, usa repr corto.
    - Para dicts/listas, usa yaml.dump “bonito” en una sola línea si es corto,
      o multilínea si hace falta.
    """
    # Escalares simples
    if isinstance(v, (int, float, str, bool)) or v is None:
        # Para strings largas, no poner comillas gigantes
        if isinstance(v, str) and len(v) > 80:
            return f'"{v[:77]}..."'
        return repr(v)

    # Estructuras: intentar una línea si cabe
    try:
        dumped = yaml.safe_dump(v, sort_keys=False, default_flow_style=True)
        dumped = dumped.strip()
        if len(dumped) <= 120 and "\n" not in dumped:
            return dumped
        # Si es largo, usar bloque multilínea legible
        dumped_block = yaml.safe_dump(
            v,
            sort_keys=False,
            default_flow_style=False,  # bloque
            allow_unicode=True,
            width=120,
            indent=2,
        ).rstrip()
        # sangrar líneas para que quede bonito en el warning
        indented = "\n    " + "\n    ".join(dumped_block.splitlines())
        return indented
    except Exception:
        return repr(v)


class ConfigLoader:
    """
    Carga un YAML, aplica valores por defecto EN MEMORIA (sin escribir el archivo)
    y emite un warning con las claves faltantes + el valor por defecto que se ha cargado.
    """

    def __init__(self, filepath="config.yaml", defaults=None, warn_missing=True):
        """
        Parameters
        ----------
        filepath : str
            Ruta al archivo YAML.
        defaults : dict | None
            Diccionario con valores por defecto esperados (estructura completa).
        warn_missing : bool
            Si True, imprime/loguea un warning con las claves completadas por defecto.
        """
        self.filepath = Path(filepath)
        self.defaults = defaults or {}
        self.warn_missing = warn_missing

        # 1) Cargar YAML (seguro ante fichero vacío o no-dict)
        if not self.filepath.exists():
            raise FileNotFoundError(f"Config file not found: {self.filepath}")

        with self.filepath.open("r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)

        if cfg is None:
            cfg = {}
        if not isinstance(cfg, dict):
            raise ValueError(f"Config file must contain a YAML mapping at top-level. Got: {type(cfg)}")

        self.config = cfg

        # 2) Aplicar defaults EN MEMORIA (no se escribe a disco)
        #    Ahora devolvemos lista de tuplas (full_key, default_value_usado)
        missing = self._apply_defaults(self.config, self.defaults)

        # 3) Advertir (sin preguntar y sin escribir el YAML)
        if self.warn_missing and missing:
            lines = ["[ConfigLoader] Faltaban claves en el YAML; se han cargado por defecto en memoria:"]
            for full_key, def_val in missing:
                lines.append(f"  - {full_key}: {_pretty_value(def_val)}")
            logging.warning("\n".join(lines))

    # ------------------------------------------------------------------
    def _apply_defaults(self, node: dict, defaults: dict, prefix: str = "") -> list[tuple[str, Any]]:
        """
        Inserta valores por defecto en 'node' si no existen (recursivo).
        Devuelve lista de tuplas (clave_completa, valor_por_defecto_insertado).
        """
        missing_pairs: list[tuple[str, Any]] = []

        for key, def_value in defaults.items():
            full_key = f"{prefix}{key}"

            if key not in node:
                node[key] = copy.deepcopy(def_value)
                missing_pairs.append((full_key, def_value))
            else:
                # Si el default es dict y el actual también: recursivo
                if isinstance(def_value, dict) and isinstance(node[key], dict):
                    missing_pairs.extend(
                        self._apply_defaults(node[key], def_value, prefix=full_key + ".")
                    )
                # Si def_value es dict pero node[key] no lo es,
                # preferimos NO reemplazar (respetar lo que viene del YAML).
                # Aquí podrías añadir validación de tipos si te interesa.
        return missing_pairs

    # ------------------------------------------------------------------
    def get(self, key, default=None):
        """Acceso de primer nivel con .get()."""
        return self.config.get(key, default)

    def __getattr__(self, name):
        """
        Acceso flexible:
          - Primero busca en primer nivel
          - Luego en diccionarios anidados a un nivel
        """
        if name in self.config:
            return self.config[name]
        for key, value in self.config.items():
            if isinstance(value, dict) and name in value:
                return value[name]
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )
    
class ConfigLoader_old:
    """ 
    Clase para cargar y acceder a parámetros de configuración desde un archivo YAML.

    Uso:
        cfg = ConfigLoader("config.yaml")
        print(cfg.filter)  # Accede al parámetro 'filter'
        print(cfg.method)  # Accede al parámetro 'method' dentro de 'alignment'
    """

    def __init__(self, filepath="config.yaml"):
        """
        Inicializa la clase cargando el archivo YAML.

        Parámetros:
            filepath (str): Ruta al archivo YAML de configuración.
        """
        with open(filepath, "r") as file:
            self.config = yaml.safe_load(file)

    def get(self, key, default=None):
        """
        Accede a una clave de primer nivel del YAML.

        Parámetros:
            key (str): Clave a buscar.
            default: Valor por defecto si no se encuentra la clave.

        Retorna:
            Valor asociado a la clave o el valor por defecto.
        """
        return self.config.get(key, default)

    def __getattr__(self, name):
        """
        Permite acceder a cualquier parámetro del YAML como atributo.

        Si el parámetro está anidado (por ejemplo, dentro de 'alignment'),
        también se puede acceder directamente por su nombre.

        Ejemplo:
            cfg.method → accede a config['alignment']['method']
        """
        if name in self.config:
            return self.config[name]
        for key, value in self.config.items():
            if isinstance(value, dict) and name in value:
                return value[name]
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")

def parse_range(value, max_value=None):
    """
    Converts strings like '1:6', '3', ':' into a list of integers.
    If max_value is provided, ':' means range(1, max_value + 1)
    """
    if isinstance(value, list):
        return value
    if isinstance(value, int):
        return [value]
    if not isinstance(value, str):
        raise ValueError(f"Unsupported type for process_ocs: {type(value)}")

    value = value.strip()

    # ':' means full range
    if value == ":":
        if max_value is None:
            raise ValueError("':' requires max_value to be defined")
        return list(range(0, max_value+1))

    # '1:6' means range 1..6
    if ":" in value:
        start, end = value.split(":")
        return list(range(int(start), int(end) + 1))

    # single number
    return [int(value)]

def timeit(method):
    """ Decorator to measure execution time of methods. """
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            t = (te - ts) * 1000

            if t > 1e3:
                logging.info(f"{method.__name__} executed in {t / 1000. :.2f} seconds")
            elif t > 1e6:
                logging.info(f"{method.__name__} executed in {t / 1000. / 60.:.2f} minutes")
            else:
                logging.info(f"{method.__name__} executed in {t :.2f} ms")

        return result

    return timed

def plt_darks(dark):
    plt.imshow(dark[0,:,:], cmap='gray',clim=(128,133))
    plt.colorbar()
    plt.title('Master Dark')
    plt.show()
    plt.imshow(dark[1,:,:], cmap='gray',clim=(128,133))
    plt.colorbar()
    plt.title('Master Dark')
    plt.show()

def plt_flats(ff_data):
    
    fig, ax = plt.subplots(2, 4, figsize=(16, 6))
    fig.tight_layout()
    for ii in range(2):
        for j in range(4):
            im = ax[ii, j].imshow(
                ff_data[0, ii, j, :,:], cmap="Greys_r"
                ,clim=(0.4,1.3))
            im.set_interpolation("none")
            # Make colorbar match the size of the panel
            cbar = plt.colorbar(im, ax=ax[ii, j], fraction=0.046, pad=0.04)
    plt.show()  # Close the figure to avoid showing it
        
def plt_level(data,roi,png_folder,name,level,label='',cclim=None):

    logging.info(f'  saving png {level} in {png_folder}/pngs/{name}_{label}.png')

    if level=='1.0':
        wn, pn, _, _ = data.shape

        fig, ax = plt.subplots(wn, pn, figsize=(16, 32))
        fig.tight_layout()
        for i in range(wn):
            plr0 = np.median(data[i, 0, roi[0]:roi[1], roi[2]:roi[3]])
            for j in range(pn):
                plr = np.median(data[i, j, roi[0]:roi[1], roi[2]:roi[3]])
                limit = (plr*0.7,plr*1.3)
                if j != 0:
                    limit = (-0.005,0.005)
                im = ax[i, j].imshow(
                    data[i, j, roi[0]:roi[1], roi[2]:roi[3]], cmap="Greys_r"
                    ,clim=limit)
                im.set_interpolation("none")
                plt.colorbar(im, ax=ax[i, j])
        plt.savefig(f"{png_folder}/pngs/{name}_{label}.png", dpi=150)
        plt.close()  # Close the figure to avoid showing it

    if level=='0.7':
        wn, pn, _, _ = data.shape

        fig, ax = plt.subplots(wn, pn, figsize=(16, 32))
        fig.tight_layout()
        for i in range(wn):
            plr0 = np.median(data[i, 0, roi[0]:roi[1], roi[2]:roi[3]])
            for j in range(pn):
                plr = np.median(data[i, j, roi[0]:roi[1], roi[2]:roi[3]])
                limit = (plr*0.7,plr*1.3)
                if j != 0:
                    limit = (-plr0*0.02+plr,plr0*0.02+plr)
                im = ax[i, j].imshow(
                    data[i, j, roi[0]:roi[1], roi[2]:roi[3]], cmap="Greys_r"
                    ,clim=limit)
                im.set_interpolation("none")
                plt.colorbar(im, ax=ax[i, j])
        plt.savefig(f"{png_folder}/pngs/{name}_{label}.png", dpi=150)
        plt.close()  # Close the figure to avoid showing it

    elif level=='0.5':

        cn, wn, pn, xs, ys = data.shape

        fig, ax = plt.subplots(wn, pn, figsize=(16, 32))
        fig.tight_layout()
        for i in range(wn):
            for j in range(pn):
                plr = np.median(data[0, i, j, roi[0]:roi[1], roi[2]:roi[3]])
                if cclim:
                    im = ax[i, j].imshow(
                    data[0, i, j, roi[0]:roi[1], roi[2]:roi[3]], cmap="Greys_r",clim=cclim)
                    # ,clim=(plr*0.7,plr*1.3))
                else:
                    im = ax[i, j].imshow(
                        data[0, i, j, roi[0]:roi[1], roi[2]:roi[3]], cmap="Greys_r")
                        # ,clim=(plr*0.7,plr*1.3))
                im.set_interpolation("none")
                plt.colorbar(im, ax=ax[i, j])
        plt.savefig(f"{png_folder}/pngs/{name}_{label}_cam0.png", dpi=150)
        plt.close()  # Close the figure to avoid showing it

        fig, ax = plt.subplots(wn, pn, figsize=(16, 32))
        fig.tight_layout()
        for i in range(wn):
            for j in range(pn):
                plr = np.median(data[1, i, j, roi[0]:roi[1], roi[2]:roi[3]])
                if cclim:
                    im = ax[i, j].imshow(
                    data[1, i, j, roi[0]:roi[1], roi[2]:roi[3]], cmap="Greys_r",clim=cclim)
                    # ,clim=(plr*0.7,plr*1.3))
                else:
                    im = ax[i, j].imshow(
                        data[1, i, j, roi[0]:roi[1], roi[2]:roi[3]], cmap="Greys_r")
                    # ,clim=(plr*0.7,plr*1.3))
                im.set_interpolation("none")
                plt.colorbar(im, ax=ax[i, j])
        plt.savefig(f"{png_folder}/pngs/{name}_{label}_cam1.png", dpi=150)
        plt.close()  # Close the figure to avoid showing it

    else:
        pass

def print_shifts_by_cam(shifts):
    """
    Print camera row/col shifts in readable tables.
    shifts: numpy array of shape (wn, 2, 2, pn)
            dims = [wn, cam, (row/col), pn]
    """
    wn, _, _, pn = shifts.shape

    for cam in range(1,2):
        print(f"\n==================== Camera {cam} ====================")

        for rc_idx, rc_label in enumerate(["row", "col"]):
            print(f"\n{rc_label.upper()} SHIFTS:")

            # Header row
            headers = ["wn\\pn"] + [f"pn{j}" for j in range(pn)]

            # Build table rows
            table = []
            for i in range(wn):
                row = [f"wn{i}"] + [f"{shifts[i, cam, rc_idx, j]:.4f}" for j in range(pn)]
                table.append(row)

            # Compute column widths
            col_widths = [max(len(h), max(len(r[c]) for r in table)) for c, h in enumerate(headers)]

            # Print header
            header_line = "  ".join(h.ljust(col_widths[i]) for i, h in enumerate(headers))
            print(header_line)
            print("-" * len(header_line))

            # Print rows
            for row in table:
                print("  ".join(row[c].ljust(col_widths[c]) for c in range(len(headers))))

def format_dict_two_rows(d, precision=6):
    """Return a compact, 2-row string of key-value pairs for logging."""
    items = list(d.items())
    half = (len(items) + 1) // 2  # split roughly in half

    def format_pair(k, v, width=20):
        # Align key left, value formatted to given precision
        if isinstance(v, (int, float)):
            return f"{k:<{width}}: {v:.{precision}f}"
        return f"{k:<{width}}: {v}"

    row1 = "  ".join(format_pair(k, v) for k, v in items[:half])
    row2 = "  ".join(format_pair(k, v) for k, v in items[half:])

    return f"{row1}\n{row2}"

import numpy as np

def _normalize_roi(shape, roi):
    """Ajusta y ordena el ROI a los límites de la imagen."""
    X, Y = shape
    if roi is None:
        return 0, X, 0, Y
    x0, x1, y0, y1 = roi
    # clamp
    x0 = max(0, min(X, x0)); x1 = max(0, min(X, x1))
    y0 = max(0, min(Y, y0)); y1 = max(0, min(Y, y1))
    # ordenar
    if x1 < x0: x0, x1 = x1, x0
    if y1 < y0: y0, y1 = y1, y0
    return x0, x1, y0, y1

def balance(imgs_cam1, imgs_cam2, roi=None, eps=1e-12, clip_percentiles=(2, 98)):
    """
    Calcula el factor de balanceo para CAM2 respecto a CAM1 en *una sola λ*.

    Si las entradas son 2D (X,Y), usa esas imágenes directamente (sin promediar).
    Si son 3D (4,X,Y), usa la media sobre modulaciones para minimizar Q,U,V.

    Parámetros
    ----------
    imgs_cam1 : np.ndarray
        (X, Y) o (4, X, Y) de cámara 1 para una λ.
    imgs_cam2 : np.ndarray
        (X, Y) o (4, X, Y) de cámara 2 para una λ.
    roi : tuple | None
        (x0, x1, y0, y1) para estimación robusta. Si None, usa toda la imagen.
    eps : float
        Tolerancia numérica para evitar divisiones por cero.
    clip_percentiles : tuple(int,int)
        Percentiles (low, high) para recortar outliers en el ratio.

    Devuelve
    --------
    scale : float
        Factor para multiplicar CAM2: cam2_bal = cam2 * scale
    gamma : float
        Estimación bruta de cam2/cam1 en el ROI (útil para logs).
    """
    # Seleccionar representaciones 2D
    if imgs_cam1.ndim == 3:
        # (4, X, Y) -> mod-average
        I1 = imgs_cam1.mean(axis=0).astype(np.float64)
    elif imgs_cam1.ndim == 2:
        I1 = imgs_cam1.astype(np.float64)
    else:
        raise ValueError(f"imgs_cam1 must be 2D or 3D, got shape={imgs_cam1.shape}")

    if imgs_cam2.ndim == 3:
        I2 = imgs_cam2.mean(axis=0).astype(np.float64)
    elif imgs_cam2.ndim == 2:
        I2 = imgs_cam2.astype(np.float64)
    else:
        raise ValueError(f"imgs_cam2 must be 2D or 3D, got shape={imgs_cam2.shape}")

    # ROI seguro
    x0, x1, y0, y1 = _normalize_roi(I1.shape, roi)
    a = I1[x0:x1, y0:y1]
    b = I2[x0:x1, y0:y1]

    # Si ROI vacío, degradar
    if a.size == 0 or b.size == 0:
        return 1.0, 1.0

    # Enmascarar pixeles con señal baja en cam1 (evitar dividir por ~0)
    abs_a = np.abs(a)
    thr = np.percentile(abs_a, 5) if abs_a.size > 0 else 0.0
    mask = abs_a > (thr + eps)

    if not np.any(mask):
        # Si no hay señal suficiente en el ROI, no balanceamos
        return 1.0, 1.0

    # Ratio robusto
    ratio = b[mask] / (a[mask] + eps)
    if ratio.size == 0 or np.all(~np.isfinite(ratio)):
        return 1.0, 1.0

    # Clip para evitar outliers
    lo, hi = clip_percentiles
    p_lo, p_hi = np.percentile(ratio[np.isfinite(ratio)], [lo, hi])
    ratio_clip = ratio[(ratio >= p_lo) & (ratio <= p_hi)]
    if ratio_clip.size == 0:
        ratio_clip = ratio[np.isfinite(ratio)]

    gamma = float(np.median(ratio_clip)) if ratio_clip.size > 0 else 1.0

    # Queremos un factor para MULTIPLICAR CAM2 completo:
    #   cam2_bal = cam2 * scale
    # Con gamma ≈ median(I2/I1), el factor correcto es:
    scale = float(1.0 / (gamma + eps))

    return scale, gamma

def _normalize_alignment_df(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    required = [
        'timestamp',
        'angle_deg','t_y','t_x','center_y','center_x',
        'scale_x','scale_y','shear_x','shear_y',
        'fun','obs_ID','line','wave'
    ]
    for c in required:
        if c not in df.columns:
            df[c] = np.nan

    num_cols = [
        'timestamp','angle_deg','t_y','t_x','center_y','center_x',
        'scale_x','scale_y','shear_x','shear_y','fun'
    ]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors='coerce')

    df['obs_ID'] = df['obs_ID'].astype(str).fillna("")
    df['line']   = df['line'].astype(str).str.strip().fillna("")
    df['wave']   = pd.to_numeric(df.get('wave', -1), errors='coerce').fillna(-1).astype(int)

    df = df.sort_values(['obs_ID','line','wave','timestamp']).reset_index(drop=True)
    return df

def _canon_line_token(x):
    if x is None:
        return None
    try:
        if isinstance(x, float) and np.isnan(x):
            return None
    except Exception:
        pass
    s = str(x).strip()
    if s == "":
        return None
    s = s.replace(",", ".")
    try:
        d = Decimal(s)
        normalized = d.normalize()
        s_fixed = format(normalized, 'f')
        if '.' in s_fixed:
            s_fixed = s_fixed.rstrip('0').rstrip('.')
        if s_fixed == "":
            s_fixed = "0"
        return s_fixed
    except InvalidOperation:
        return s

def interpolate_filter(
    df: pd.DataFrame,
    timestamp: float,
    line,
    obs_ID,
    wave: int | None = None,
    default_values: dict | None = None,
    *,
    max_dt: float | None = None  # tolerancia máxima (mismas unidades que 'timestamp'); None -> sin límite
):
    """
    Selección por 'timestamp' MÁS CERCANO (no interpola):
      - Filtra por (obs_ID, line canónica, wave).
      - Si wave es int: busca en esa wave; si no hay, fallback a wave=-1.
      - Si wave es None (dataset): prioriza wave=-1; si no hay, usa todas las filas filtradas.
      - Elige SIEMPRE la fila cuyo 'timestamp' es más cercano al solicitado y registra la wave efectiva usada.
    """
    if default_values is None:
        default_values = {
            'angle_deg': 0.0, 't_y': 0.0, 't_x': 0.0,
            'center_y': 0.0, 'center_x': 0.0,
            'scale_x': 1.0, 'scale_y': 1.0,
            'shear_x': 0.0, 'shear_y': 0.0
        }

    # 1) Normalizar DF
    df = _normalize_alignment_df(df)

    # 2) Canonizar 'line'
    if 'line_key' not in df.columns:
        df['line_key'] = df['line'].apply(_canon_line_token)
    line_key = _canon_line_token(line)

    # 3) Filtro base por obs_ID + line
    df_f = df[
        (df['obs_ID'].astype(str).str.strip() == str(obs_ID).strip()) &
        (df['line_key'] == line_key)
    ].copy()

    if df_f.empty:
        logging.warning(
            f"[interpolate_filter] Sin filas para obs_ID='{obs_ID}', line='{line}' (canon='{line_key}'). Devuelvo defaults."
        )
        return dict(default_values)

    # 4) Filtro por wave + fallback, y determinar 'used_wave' que vamos a reportar
    def _by_wave(dfi, w):
        return dfi[dfi['wave'].astype(int) == int(w)]

    df_w = None
    used_wave: int | str

    if wave is not None:
        df_w = _by_wave(df_f, wave)
        if df_w.empty:
            logging.info(f"[interpolate_filter] Sin datos para wave={wave}. Fallback a wave=-1 (global).")
            df_w = _by_wave(df_f, -1)
            used_wave = -1 if not df_w.empty else "mixed"  # si tampoco hay -1, luego usaremos todas las filas base
        else:
            used_wave = int(wave)
    else:
        # Dataset-level: prioriza wave == -1
        df_w = _by_wave(df_f, -1)
        if df_w.empty:
            logging.info("[interpolate_filter] Dataset-level sin wave=-1. Uso todas las filas base (mezcla de waves).")
            df_w = df_f
            used_wave = "mixed"  # se resolverá con la fila elegida
        else:
            used_wave = -1

    if df_w.empty:
        logging.warning(f"[interpolate_filter] Sin datos tras aplicar wave (solicitado={wave}). Devuelvo defaults.")
        return dict(default_values)

    # 5) Elegir SIEMPRE la fila con timestamp MÁS CERCANO
    df_w = df_w.dropna(subset=['timestamp']).sort_values('timestamp')
    t_vals = df_w['timestamp'].values.astype(float)
    idx = int(np.argmin(np.abs(t_vals - float(timestamp))))
    row = df_w.iloc[idx]
    dt = abs(float(row['timestamp']) - float(timestamp))

    # Si 'used_wave' era "mixed" o no teníamos claro el fallback, tomar la wave real de la fila elegida
    # Esto cubre:
    #   - dataset-level sin -1 (se usaron todas y ahora conocemos la wave elegida)
    #   - fallback a -1 inexistente (si hubiéramos caído a base y mezclado)
    if used_wave == "mixed":
        try:
            used_wave = int(row.get('wave', -1))
        except Exception:
            used_wave = -1

    logging.info(
        "[interpolate_filter] elegido timestamp=%.6f (Δt=%.6f) para obs_ID=%s line=%s wave_usada=%s",
        float(row['timestamp']), dt, obs_ID, line, str(used_wave)
    )

    if (max_dt is not None) and (dt > max_dt):
        logging.warning(
            f"[interpolate_filter] Vecino más cercano a t={timestamp:.6f} está demasiado lejos (Δt={dt:.6f} > max_dt={max_dt}). Devuelvo defaults."
        )
        return dict(default_values)

    # 6) Devolver la fila elegida (sin interpolar)
    return {
        'angle_deg': float(row['angle_deg']),
        't_y':       float(row['t_y']),
        't_x':       float(row['t_x']),
        'center_y':  float(row['center_y']),
        'center_x':  float(row['center_x']),
        'scale_x':   float(row['scale_x']),
        'scale_y':   float(row['scale_y']),
        'shear_x':   float(row['shear_x']),
        'shear_y':   float(row['shear_y']),
    }

def parse_header_time(val):
    if isinstance(val, datetime):
        return val
    if isinstance(val, str):
        for fmt in ("%d%m%YT%H%M%S.%f", "%d%m%YT%H%M%S"):
                try:
                    return datetime.strptime(val, fmt)
                except Exception:
                    pass
    try:
        return datetime.fromtimestamp(float(val))
    except Exception:
        raise ValueError(f"Unsupported header time format: {val!r}")

def minutes_from_dt(d):
    return d.day * 1440 + d.hour * 60 + d.minute + d.second/60.0 + d.microsecond/60000000.0

def _timestamp_from_filename(filename):
    # Busca patrón tipo 10072024T191616
    match = re.search(r'(\d{8}T\d{6})', filename)
    if not match:
        raise ValueError("No date-time pattern found in filename.")
    
    datetime_str = match.group(1)
    date_obj = datetime.strptime(datetime_str, "%d%m%YT%H%M%S")

    # Timestamp en minutos desde inicio de mes (tu definición)
    timestamp = date_obj.day * 1440 + date_obj.hour * 60 + date_obj.minute  

    return int(timestamp)


def remove_freq(image, kx, ky, h):
    ny, nx = image.shape
    cy, cx = 800, 800

    temp = np.full((1600, 1600), np.median(image))
    temp[:ny, :nx] = image

    fft = np.fft.fftshift(np.fft.fft2(temp))

    for kx_, ky_ in zip(kx, ky):
        fft[ky_ - h + cy : ky_ + h + 1 + cy,
            kx_ - h + cx : kx_ + h + 1 + cx] = 0

        fft[-ky_ - h + cy : -ky_ + h + 1 + cy,
            -kx_ - h + cx : -kx_ + h + 1 + cx] = 0

    fft = np.fft.ifft2(np.fft.ifftshift(fft))

    return np.real(fft[:ny, :nx])




def correct_image_fft(
    img,
    kx,
    ky,
    h=3,
    thr=0.003,
    mode="stokes",
    sigma_bg=50,
):
    """
    FFT fringe correction for a single 2D image.

    mode = "stokes"    -> assumes signal centered at zero, uses threshold
    mode = "intensity" -> removes background, filters full residual (no threshold)
    """

    if mode == "intensity":
        # remove slow background (offset)
        offset = gaussian_filter(img, sigma=sigma_bg)
        residual = img - offset

        # filter EVERYTHING in residual
        filtered = remove_freq(residual, kx, ky, h)

        # restore background
        return filtered + offset

    else:  # "stokes"
        low  = np.where(np.abs(img) > thr, 0, img)
        high = img - low

        return high + remove_freq(low, kx, ky, h)

def correct_data_fft(
    data,
    kx,
    ky,
    h=3,
    thr=0.003,
    mode="stokes",
    sigma_bg=50,
):
    """
    Apply FFT fringe correction to N-D data.
    Assumes last two axes are (y, x).
    """

    data_ = data.copy()
    leading_shape = data_.shape[:-2]

    for idx in np.ndindex(leading_shape):
        img = data_[idx + (slice(None), slice(None))]

        # # Calculate power spectrum with apodization and the inverse Fourier transform
        # f_fft, power_spectrum = calculate_power_spectrum_with_apodization(data_[idx + (slice(None), slice(None))], window_type='hanning')

        # # go back
        # inverse_image = calculate_inverse_fourier_transform(f_fft)

        # # Plot the results principales
        # fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        # axs[0].imshow(data_[idx + (slice(None), slice(None))], cmap='gray')
        # axs[0].set_title('Original Image')
        # axs[1].imshow(np.log10(power_spectrum[800-40:800+40,800-40:800+40]), cmap='gray')
        # axs[1].set_title('Power Spectrum')
        # axs[2].imshow(inverse_image.real, cmap='gray')
        # axs[2].set_title('Inverse Fourier Transform')
        # plt.show()
        
        data_[idx + (slice(None), slice(None))] = correct_image_fft(
            img,
            kx,
            ky,
            h=h,
            thr=thr,
            mode=mode,
            sigma_bg=sigma_bg,
        )

        # # Calculate power spectrum with apodization and the inverse Fourier transform
        # f_fft, power_spectrum = calculate_power_spectrum_with_apodization(data_[idx + (slice(None), slice(None))], window_type='hanning')

        # # go back
        # inverse_image = calculate_inverse_fourier_transform(f_fft)

        # # Plot the results principales
        # fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        # axs[0].imshow(data_[idx + (slice(None), slice(None))], cmap='gray')
        # axs[0].set_title('Original Image')
        # axs[1].imshow(np.log10(power_spectrum[800-40:800+40,800-40:800+40]), cmap='gray')
        # axs[1].set_title('Power Spectrum')
        # axs[2].imshow(inverse_image.real, cmap='gray')
        # axs[2].set_title('Inverse Fourier Transform')
        # plt.show()


    return data_

def calculate_power_spectrum_with_apodization(image, window_type=None):
    """
    Calculate the power spectrum of an image with apodization to avoid artifacts.

    Parameters:
        image (ndarray): Input 2D image array.
        window_type (str): Type of window to apply for apodization. Default is 'hanning'.
                          Other options include 'hamming', 'blackman', etc.

    Returns:
        power_spectrum (ndarray): 2D power spectrum of the image after applying apodization.
    """
    
    # Get the size of the image
    nx, ny = image.shape
    
    # Create a 1D window (Hanning by default, can be changed)
    if window_type:
        if window_type == 'hanning':
            window_x = np.hanning(nx)
            window_y = np.hanning(ny)
        elif window_type == 'hamming':
            window_x = np.hamming(nx)
            window_y = np.hamming(ny)
        elif window_type == 'blackman':
            window_x = np.blackman(nx)
            window_y = np.blackman(ny)
        else:
            raise ValueError(f"Unsupported window type: {window_type}")
        
        # Create the 2D window by outer product
        window_2d = np.outer(window_x, window_y)
        
        # Apply the window to the image (element-wise multiplication)
        apodized_image = image * window_2d
    else:
        apodized_image = np.copy(image)

    # Perform 2D FFT (shift before and after for proper centering)
    f_fft = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(apodized_image)))
    
    # Compute the power spectrum (magnitude squared of the FFT)
    power_spectrum = np.abs(f_fft) ** 2
    
    return f_fft, power_spectrum

def calculate_inverse_fourier_transform(f_fft):
    """
    Calculate the inverse Fourier transform of a given 2D FFT.

    Parameters:
        f_fft (ndarray): Input 2D FFT array.

    Returns:
        inverse_image (ndarray): 2D image obtained from the inverse Fourier transform.
    """
    
    # Perform inverse 2D FFT (shift before and after for proper centering)
    inverse_image = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(f_fft)))
    
    # Take the real part of the inverse image
    inverse_image = np.real(inverse_image)
    
    return inverse_image


def generar_cuadrantes_nx_n(H, W, n):
    """
    Divide una imagen en una cuadrícula n×n de cuadrados del mismo tamaño.

    Parámetros:
        H, W : int
            Alto (H) y ancho (W) de la imagen.
        n : int
            Número de cuadrados por eje (n=2 -> 2x2, n=3 -> 3x3, etc.)

    Retorna:
        Lista de tuplas (y1, y2, x1, x2) que representan las coordenadas
        de cada cuadrado en formato (fila_superior, fila_inferior, col_izquierda, col_derecha).
    """
    # tamaño ideal de cada cuadrado
    base_h = H / n
    base_w = W / n

    rois = []
    for i in range(n):
        for j in range(n):
            y1 = int(round(i * base_h))
            y2 = int(round((i + 1) * base_h))
            x1 = int(round(j * base_w))
            x2 = int(round((j + 1) * base_w))
            rois.append((y1, y2, x1, x2))

    return rois


def extract_crosstalk_coeffs(cfg_crosstalk: dict):
    """
    Extrae de cfg['crosstalk']:
      - flags de aplicación (apply_only, use_local, local_order, deriv_sigma, channels)
      - coeficientes: global o per_wavelength
      - per_wavelength: bool

    Devuelve:
      params : dict {apply_only, use_local, local_order, deriv_sigma, channels}
      coeffs : dict global {'Q':{a,b,c,d,e},...} o lista por λ [ {'Q':{...},...}, ... ]
      per_wavelength : bool
    """

    def _coerce_float(x, default=0.0):
        try:
            return float(x)
        except Exception:
            return float(default)

    def _normalize_channels(chs):
        # Acepta lista/tupla de strings y filtra a ('Q','U','V') en ese orden
        valid = ('Q','U','V')
        if chs is None:
            return valid
        out = []
        for c in chs:
            c2 = str(c).strip().upper()
            if c2 in valid:
                out.append(c2)
        return tuple(out) if out else valid

    def _ensure_block_coeffs(block, channels=('Q','U','V')):
        """
        Asegura que 'block[ch]' existe y tiene las claves a,b,c,d,e numéricas.
        Modifica 'block' in-place.
        """
        for ch in channels:
            if ch not in block or block[ch] is None:
                block[ch] = {}
            for k in ('a','b','c','d','e'):
                block[ch][k] = _coerce_float(block[ch].get(k, 0.0), 0.0)

    if not isinstance(cfg_crosstalk, dict):
        raise ValueError("cfg['crosstalk'] debe ser un diccionario.")

    # ---- parámetros/flags ----
    apply_only  = bool(cfg_crosstalk.get('apply_only', False))
    use_local   = bool(cfg_crosstalk.get('use_local', False))
    local_order = int(cfg_crosstalk.get('local_order', 1))
    deriv_sigma = float(cfg_crosstalk.get('deriv_sigma', 0.0))
    aggregate_wavelengths = bool(cfg_crosstalk.get('aggregate_wavelengths', False))
    channels    = _normalize_channels(cfg_crosstalk.get('channels', ('Q','U','V')))

    # ---- coeficientes ----
    coeffs_global = cfg_crosstalk.get('coeffs_global', None)
    coeffs_per_wv = cfg_crosstalk.get('coeffs_per_wavelength', None)

    if coeffs_global is None and coeffs_per_wv is None:
        raise ValueError(
            "En cfg['crosstalk'] falta 'coeffs_global' o 'coeffs_per_wavelength'."
        )

    # Forma de los coeficientes
    per_wavelength = coeffs_per_wv is not None
    coeffs = coeffs_per_wv if per_wavelength else coeffs_global

    # ---- validación / coerción numérica ----
    if not per_wavelength:
        if not isinstance(coeffs, dict):
            raise ValueError("'coeffs_global' debe ser dict por canal.")
        _ensure_block_coeffs(coeffs, channels=channels)
    else:
        if not isinstance(coeffs, (list, tuple)):
            raise ValueError("'coeffs_per_wavelength' debe ser una lista/tupla de dicts por λ.")
        for i, block in enumerate(coeffs):
            if not isinstance(block, dict):
                raise ValueError(f"Elemento #{i} en 'coeffs_per_wavelength' no es dict.")
            _ensure_block_coeffs(block, channels=channels)

    params = dict(
        apply_only=apply_only,
        use_local=use_local,
        local_order=local_order,
        deriv_sigma=deriv_sigma,
        aggregate_wavelengths=aggregate_wavelengths,
        channels=channels
    )
    return params, coeffs, per_wavelength
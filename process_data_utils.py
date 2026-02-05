import yaml
import logging
import time
import matplotlib.pyplot as plt
import numpy as np

class ConfigLoader:
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
        
def plt_level(data,roi,png_folder,name,level,label=''):

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

    for cam in range(2):
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

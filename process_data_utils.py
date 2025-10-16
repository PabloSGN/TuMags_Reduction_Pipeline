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

    if level=='1.0':
        wn, pn, _, _ = data.shape

        fig, ax = plt.subplots(wn, pn, figsize=(16, 32))
        fig.tight_layout()
        for i in range(wn):
            for j in range(pn):
                plr = np.median(data[i, j, roi[0]:roi[1], roi[2]:roi[3]])
                if j == 0:
                    levels = (0.2,1.2)
                else:
                    levels = (0.01,0.01)

                im = ax[i, j].imshow(
                    data[i, j, roi[0]:roi[1], roi[2]:roi[3]], cmap="Greys_r"
                    ,clim=levels)
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
                    data[0, i, j, roi[0]:roi[1], roi[2]:roi[3]], cmap="Greys_r"
                    ,clim=(plr*0.7,plr*1.3))
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
                    data[1, i, j, roi[0]:roi[1], roi[2]:roi[3]], cmap="Greys_r"
                    ,clim=(plr*0.7,plr*1.3))
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
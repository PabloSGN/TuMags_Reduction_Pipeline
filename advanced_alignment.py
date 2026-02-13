
from scipy.ndimage import affine_transform
from scipy.optimize import minimize
import matplotlib.patches as patches
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import logging
from process_data_utils import balance
from demodulation import demodulate
import time

# import (ConfigLoader,parse_range, 
#                               print_shifts_by_cam, plt_darks, 
#                               plt_flats,plt_level,
#                               balance,
#                               interpolate_filter) #, format_dict_two_rows

def _generate_sampling_points(image_shape, margin=200):
    """
    Devuelve 9 puntos:
    0–3 : esquinas
    4–7 : centros de los lados
    8   : centro
    Formato (y, x).
    """
    H, W = image_shape[:2]
    cy, cx = H // 2, W // 2

    points = [
        # Esquinas
        (margin, margin),            # 0 top-left
        (margin, W - margin),        # 1 top-right
        (H - margin, margin),        # 2 bottom-left
        (H - margin, W - margin),    # 3 bottom-right

        # Centros de los lados
        (margin, cx),                # 4 top-center
        (H - margin, cx),            # 5 bottom-center
        (cy, margin),                # 6 left-center
        (cy, W - margin),            # 7 right-center

        # Centro
        (cy, cx)                     # 8 center
    ]
    return points


def apply_transform(image, angle_deg, t, center,
                    scale_x=1.0, scale_y=1.0,
                    shear_x=0.0, shear_y=0.0):

    angle_rad = np.radians(angle_deg)
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)

    R = np.array([
        [cos_a, -sin_a],
        [sin_a,  cos_a]
    ])

    Distortion = np.array([
        [scale_x, shear_x],
        [shear_y, scale_y]
    ])

    A = Distortion @ R
    A_inv = np.linalg.inv(A)

    offset = center - A_inv @ (center + t)

    return affine_transform(
        image,
        matrix=A_inv,
        offset=offset,
        order=3,
        mode='nearest',
        cval=0.0
    )

def _clip_center_to_image(y, x, half, H, W):
    """
    Asegura que el cuadrado de lado 2*half centrado en (y,x) queda dentro de [0..H) x [0..W).
    Si ya cabe, no modifica nada.
    """
    y = int(np.clip(y, half, H - 1 - half))
    x = int(np.clip(x, half, W - 1 - half))
    return y, x

def _extract_patch(image, center, size=20):
    """
    Extrae un parche cuadrado centrado en `center` (y,x) con lado `size`, 
    forzando a que quede dentro de la imagen.
    """
    y, x = np.round(center).astype(int)
    half = size // 2
    H, W = image.shape[:2]

    # Clip de seguridad. Si half <= margin, no moverá nada.
    y, x = _clip_center_to_image(y, x, half, H, W)

    return image[y - half : y + half, x - half : x + half]

_call_counter = 0

def _cost_function_pixelwise(params, image0, image1,
                            size_corner=300,
                            size_center=180, _call_counter = _call_counter):

    # global _call_counter
    _call_counter += 1

    angle_deg, t_y, t_x, c_y, c_x, scale_x, scale_y, shear_x, shear_y = params

    t = np.array([t_y, t_x])
    center = np.array([c_y, c_x])

    transformed_image1 = apply_transform(
        image1, angle_deg, t, center,
        scale_x, scale_y, shear_x, shear_y
    )

    points = _generate_sampling_points(image0.shape)

    weighted_rms = []
    weights = []

    for i, p in enumerate(points):

        # --- tamaños ---
        if i == 8:        # centro
            size = size_center
            weight = 0.5
        elif i < 4:       # esquinas
            size = size_corner
            weight = 2.0
        else:             # lados
            size = size_corner
            weight = 1.0

        patch0 = _extract_patch(image0, p, size=size)
        patch1 = _extract_patch(transformed_image1, p, size=size)

        if patch0.shape == patch1.shape:
            diff = patch0 - patch1
            rms = np.sqrt(np.mean(diff**2))
            weighted_rms.append(weight * rms)
            weights.append(weight)

    mean_rms = np.sum(weighted_rms) / np.sum(weights)

    # ---- DEBUG PRINT ----
    # if _call_counter % 10 == 0:
    print(
        f"[{_call_counter:04d}] "
        f"angle={angle_deg:+.6f}  "
        f"ty={t_y:+.2f} tx={t_x:+.2f}  "
        f"RMS={mean_rms:.6f}  "
        f"cy={c_y:+.2f} cx={c_x:+.2f}  "

    )

    return mean_rms

def plot_sampling_patches(image1, image2, size_corner=100, size_center=180, margin=200):
    """
    Dibuja los 9 parches en `image1` garantizando que queden dentro.
    Los tamaños y el margen se validan para que el dibujo tenga sentido.
    """
    H, W = image1.shape[:2]

    # --- Validaciones coherentes con la regla "el cuadrado dentro de la imagen" ---
    half_corner = size_corner // 2
    half_center = size_center // 2

    # Para esquinas y lados (que están a 'margin' del borde), necesitamos:
    if half_corner > margin:
        raise ValueError(
            f"size_corner={size_corner} es demasiado grande para margin={margin} "
            f"(necesitas size_corner/2 <= margin)."
        )

    # Comprobación básica de que los parches caben también por el otro lado
    if any([
        margin - half_corner < 0,               # arriba/izquierda
        W - margin - half_corner < 0,           # derecha
        H - margin - half_corner < 0            # abajo
    ]):
        raise ValueError("El tamaño del parche o el margen hace que no quepan dentro de la imagen.")

    points = _generate_sampling_points(image1.shape, margin)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(image1-image2, cmap='gray')
    ax.set_title("Sampling patches (weights encoded)")
    ax.axis('off')

    for i, (y, x) in enumerate(points):
        # --- Determinar tamaño/color/etiqueta ---
        if i == 8:          # centro
            size = size_center
            color = 'red'
            label = 'center'
        elif i < 4:         # esquinas
            size = size_corner
            color = 'yellow'
            label = 'corner'
        else:               # lados
            size = size_corner
            color = 'lime'
            label = 'side'

        half = size // 2

        # Clip de seguridad (no moverá si half <= margin)
        y_safe, x_safe = _clip_center_to_image(y, x, half, H, W)

        # Dibujo del rectángulo
        rect = patches.Rectangle(
            (x_safe - half, y_safe - half),
            size, size,
            linewidth=2,
            edgecolor=color,
            facecolor='none'
        )
        ax.add_patch(rect)

        # Índice del punto
        ax.text(
            x_safe, y_safe, f"{i}",
            color=color,
            fontsize=10,
            ha='center',
            va='center'
        )

        # (Opcional) pequeña cruz en el centro real solicitado y otro en el centro ajustado:
        ax.plot([x], [y], marker='+', color=color, markersize=8, mew=2, alpha=0.6)        # pedido
        if (y_safe != y) or (x_safe != x):
            ax.plot([x_safe], [y_safe], marker='x', color=color, markersize=6, mew=2)     # ajustado

    plt.show()

def _print_result(res):
    names = [
        "angle_deg",
        "t_y",
        "t_x",
        "center_y",
        "center_x",
        "scale_x",
        "scale_y",
        "shear_x",
        "shear_y"
    ]

    print("\n=== RESULTADO FINAL DE OPTIMIZACIÓN ===")
    print(f"Convergió: {res.success}")
    print(f"Mensaje:   {res.message}")
    print("\nParámetros óptimos:\n")

    for name, value in zip(names, res.x):
        print(f"  {name:10s}: {value: .6f}")

    print(f"\nCosto final: {res.fun:.6f}")
    print("=======================================\n")

def align_advance(I_cam1, I_cam2):

    size_corner = 500   # half=150
    size_center = 300   # half=100
    margin = 200
    # plot_sampling_patches(
    #     I_cam1, I_cam2,
    #     size_corner=size_corner,
    #     size_center=size_center,
    #     margin=margin
    # )
    args=(I_cam1, I_cam2,size_corner,size_center,_call_counter)

    init_params = [
        0.08,   # angle_deg
        0.1,    # t_y
        0.1,    # t_x
        800,    # center_y
        800,    # center_x
        1.0,    # scale_x
        1.0,    # scale_y
        0.0,    # shear_x
        0.0     # shear_y
    ]


    H, W = I_cam1.shape[:2]

    bounds = [
        (-1.0, 1.0),        # angle_deg
        (-5.0,  5.0),       # t_y
        (-5.0,  5.0),       # t_x
        (200,     1400),         # center_y
        (200,     1400),         # center_x
        (0.99, 1.01),     # scale_x (FIJO)
        (0.99, 1.01),     # scale_y (FIJO)
        (-0.01, 0.01),     # shear_x (FIJO)
        (-0.01, 0.01),     # shear_y (FIJO)
    ]

    res = minimize(
        _cost_function_pixelwise,
        init_params,
        args=args,
        method='Powell',
        bounds=bounds,
        options={
            'maxiter': 100,
            'disp': True,
        }
    )

    _print_result(res)

    return res


def update_alignment_csv(csv_path, timestamp, parameter_vector, cost, overwrite = False):

    csv_path = Path(csv_path)

    col_names = [
            "timestamp",
            "angle_deg", "t_y", "t_x",
            "center_y", "center_x",
            "scale_x", "scale_y",
            "shear_x", "shear_y",
            "cost"
    ]

    # Convert result into row dict
    row = {
            "timestamp": timestamp,
            "angle_deg": parameter_vector[0],
            "t_y": parameter_vector[1],
            "t_x": parameter_vector[2],
            "center_y": parameter_vector[3],
            "center_x": parameter_vector[4],
            "scale_x": parameter_vector[5],
            "scale_y": parameter_vector[6],
            "shear_x": parameter_vector[7],
            "shear_y": parameter_vector[8],
            "cost": cost
    }

    # Si no existe → crear
    if not csv_path.exists():
            df = pd.DataFrame([row], columns=col_names)
            df.to_csv(csv_path, index=False)
            return

    # Si sí existe → cargar y reemplazar
    df = pd.read_csv(csv_path)

    # Buscar si el timestamp ya existe
    if timestamp in df["timestamp"].values and not overwrite:
            df.loc[df["timestamp"] == timestamp] = list(row.values())
    else:
            df.loc[len(df)] = list(row.values())

    # Guardar actualizado
    df.to_csv(csv_path, index=False)

def advance_alignment(data, filter, wave):
    
    tic = time.time()

    _, demod = demodulate(data[:, wave], 
                            filt=filter, 
                            onelambda=True, 
                            BothCams=True)

    I_cam1 = demod[0, 0]
    I_cam2 = demod[1, 0]

    scale, gamma = balance(I_cam1, I_cam2)
    I_cam2 = I_cam2 * scale
    logging.info(f"Balance (in advance_alignment): wave: {wave}, scale={scale:.6f}, gamma={gamma:.6f}")

    # logging.info("Global alignment between cameras (in advanced mode) using Stokes I")
    # a,b = np.gradient(I_cam1)
    # I_cam1 = a+b
    # a,b = np.gradient(I_cam2)
    # I_cam2 = a+b

    result = align_advance(I_cam1,I_cam2)
    angle_deg, t_y, t_x, c_y, c_x, scale_x, scale_y, shear_x, shear_y = result.x
    t = np.array([t_y, t_x])
    center = np.array([c_y, c_x])

    # Aplicar shift global SOLO a cam2 
    for ld in range(data.shape[1]):
        for j in range(data.shape[2]):
            data[1, ld, j] = apply_transform(
                data[1, ld, j], angle_deg, t, center,
                scale_x, scale_y, shear_x, shear_y
            )

    tac = time.time()

    logging.info(f"Alignment finished in {round(tac - tic, 3)} s.")

    return data, result

# ---------------------------- DESCRIPTION --------------------------------------- #
"""

Module with all the functions related to the alignment of the observation modes. 

Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ IMPORTS ----------------------------------------- #

import numpy as np
import time
import logging
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.fftpack import fftshift, ifftshift, fft2, ifft2

# Own functions
from demodulation import demodulate
from demodulation import mod_matrices_david_ct as mod_matrices
from process_data_utils import balance

# ------------------------------  AUX FUNCTIONS  --------------------------------- # 

def plot_maps(dy, dx, title=""):
    vmax = max(np.abs(dy).max(), np.abs(dx).max())
    vmin = -vmax

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))

    im0 = axs[0].imshow(dy, cmap='coolwarm', vmin=vmin, vmax=vmax)
    axs[0].set_title(title + " dy")
    plt.colorbar(im0, ax=axs[0])

    im1 = axs[1].imshow(dx, cmap='coolwarm', vmin=vmin, vmax=vmax)
    axs[1].set_title(title + " dx")
    plt.colorbar(im1, ax=axs[1])

    plt.tight_layout()
    plt.show()

def build_displacement_maps_simple(H, W, quadrants_roi, shifts_q):
    """
    Crea mapas dy y dx por píxel asignando el desplazamiento de cada cuadrante
    a todo su ROI, sin suavizado ni mezcla.
    """
    dy_map = np.zeros((H, W), dtype=float)
    dx_map = np.zeros((H, W), dtype=float)

    for (y1, y2, x1, x2), (dy, dx) in zip(quadrants_roi, shifts_q):
        dy_map[y1:y2, x1:x2] = dy
        dx_map[y1:y2, x1:x2] = dx

    return dy_map, dx_map

def _percentile_limits(img, plow=1, phigh=99):
    vmin = np.percentile(img, plow) if img.size else 0.0
    vmax = np.percentile(img, phigh) if img.size else 1.0
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin, vmax = float(np.nanmin(img)), float(np.nanmax(img)) if img.size else (0.0, 1.0)
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
            vmin, vmax = 0.0, 1.0
    return vmin, vmax

def _apply_grid(ax, step=10, grid_color="cyan", lw=0.5, alpha=0.7):
    """Configura una rejilla visible sin mostrar etiquetas."""
    ax.xaxis.set_major_locator(MultipleLocator(step))
    ax.yaxis.set_major_locator(MultipleLocator(step))
    ax.grid(True, which='major', color=grid_color, linewidth=lw, alpha=alpha)
    # Oculta etiquetas y marcas, pero conserva las posiciones de ticks
    ax.tick_params(axis='both', which='both',
                   labelbottom=False, labelleft=False,
                   bottom=False, left=False, length=0)

def _cosine_ramp(n):
    # rampa coseno de 0->1 en n muestras (n>=1)
    if n <= 0:
        return np.array([], dtype=np.float64)
    t = np.linspace(0, np.pi, n, endpoint=True)
    return 0.5 - 0.5*np.cos(t)  # 0..1

def _feather_mask_from_core(h, w, core, margin_y, margin_x):
    """
    Crea una máscara 2D con plateau=1 en 'core' (en coords del parche)
    y rampas coseno en los márgenes hacia fuera (hasta 'margin_*').
    Fuera de core + rampas -> 0.
    """
    y0, y1, x0, x1 = core
    y0 = max(0, min(h, y0)); y1 = max(0, min(h, y1))
    x0 = max(0, min(w, x0)); x1 = max(0, min(w, x1))

    wy = np.zeros(h, dtype=np.float64)
    wx = np.zeros(w, dtype=np.float64)

    # Eje Y
    mY_top = min(margin_y, y0)              # espacio para rampa arriba
    mY_bot = min(margin_y, h - y1)          # espacio para rampa abajo
    if y0 - mY_top < y1 + mY_bot:
        # tramo superior: 0 -> 1
        if mY_top > 0:
            wy[y0 - mY_top:y0] = _cosine_ramp(mY_top)
        # plateau
        wy[y0:y1] = 1.0
        # tramo inferior: 1 -> 0
        if mY_bot > 0:
            wy[y1:y1 + mY_bot] = _cosine_ramp(mY_bot)[::-1]

    # Eje X
    mX_left  = min(margin_x, x0)
    mX_right = min(margin_x, w - x1)
    if x0 - mX_left < x1 + mX_right:
        if mX_left > 0:
            wx[x0 - mX_left:x0] = _cosine_ramp(mX_left)
        wx[x0:x1] = 1.0
        if mX_right > 0:
            wx[x1:x1 + mX_right] = _cosine_ramp(mX_right)[::-1]

    mask = np.outer(wy, wx)
    return mask

def _cosine_taper_2d(h, w, margin_y, margin_x):
    # Hann 2D sobre los bordes del parche (pre-apodización, pequeña)
    my = int(max(0, margin_y)); mx = int(max(0, margin_x))
    wy = np.ones(h, dtype=np.float64)
    wx = np.ones(w, dtype=np.float64)
    if my > 0:
        r = _cosine_ramp(my)
        wy[:my] = r
        wy[-my:] = r[::-1]
    if mx > 0:
        r = _cosine_ramp(mx)
        wx[:mx] = r
        wx[-mx:] = r[::-1]
    return np.outer(wy, wx)

def apply_quadrant_shifts_with_feather(
    img2d,                       # imagen 2D (X,Y)
    quadrants_roi,               # lista de (y1,y2,x1,x2)
    shifts_q,                    # lista de (dy,dx) por cuadrante (en mismo orden)
    margin=10,                   # px de margen para feather (solape)
    shift_func=None,             # callable(patch, shift=[dy,dx], wrap=False, fill=0)
    wrap=False, fill=0,
    pre_apod=3                   # apodización pequeña pre-shift (px)
):
    """
    - Expande cada ROI en margin + |shift| (por eje) para garantizar cobertura tras el desplazamiento.
    - Aplica ventana coseno suave en los bordes del parche (pre_apod) antes del shift (reduce ringing).
    - Construye una máscara de escritura con plateau=1 exactamente en el ROI original y rampas en el solape.
    - Mezcla y normaliza (suma ponderada); partición de unidad en solapes -> sin “cruz” en el centro.
    """
    H, W = img2d.shape
    num = np.zeros_like(img2d, dtype=np.float64)
    den = np.zeros_like(img2d, dtype=np.float64)

    for (y1, y2, x1, x2), (dy, dx) in zip(quadrants_roi, shifts_q):
        # expansión dependiente del shift (garantiza cobertura post-shift)
        pad_y = margin + int(np.ceil(abs(dy)))
        pad_x = margin + int(np.ceil(abs(dx)))

        ye1 = max(0, y1 - pad_y)
        ye2 = min(H, y2 + pad_y)
        xe1 = max(0, x1 - pad_x)
        xe2 = min(W, x2 + pad_x)
        if ye2 <= ye1 or xe2 <= xe1:
            continue

        patch = img2d[ye1:ye2, xe1:xe2]
        h, w = patch.shape

        # 1) Apodización pre-shift (pequeña, sobre bordes del parche)
        apy = min(pre_apod, h//2 - 1) if h >= 3 else 0
        apx = min(pre_apod, w//2 - 1) if w >= 3 else 0
        pre_w = _cosine_taper_2d(h, w, apy, apx)
        patch_win = patch * pre_w

        # 2) Desplaza el parche
        shifted = shift_func(patch_win, shift=[dy, dx], wrap=wrap, fill=fill)

        # 3) Máscara de escritura basada en el núcleo (ROI original)
        core_y0 = y1 - ye1
        core_y1 = y2 - ye1
        core_x0 = x1 - xe1
        core_x1 = x2 - xe1

        # márgenes efectivos (no más grandes que el espacio disponible)
        my = min(margin, h//2)
        mx = min(margin, w//2)
        write_w = _feather_mask_from_core(
            h, w,
            core=(core_y0, core_y1, core_x0, core_x1),
            margin_y=my, margin_x=mx
        )

        # 4) Acumulación
        num[ye1:ye2, xe1:xe2] += shifted * write_w
        den[ye1:ye2, xe1:xe2] += write_w

    out = img2d.copy()
    m = den > 1e-12
    out[m] = num[m] / den[m]
    return out

def apply_quadrant_shifts_with_feather_stack(
    stack4, quadrants_roi, shifts_q_per_mod, margin=10, shift_func=None, wrap=False, fill=0
):
    """
    Aplica shifts por cuadrante a un stack (4,X,Y).
    shifts_q_per_mod: lista (por mod j) de listas (por cuadrante q) de tuples (dy,dx),
                      es decir: longitud 4, cada elemento longitud Q.
    """
    out = stack4.copy()
    nmods, H, W = stack4.shape
    assert len(shifts_q_per_mod) == nmods, "Se espera una lista de shifts por mod (len=4)."

    for j in range(nmods):
        shifts_q = shifts_q_per_mod[j]  # lista de (dy,dx) por cuadrante
        out[j] = apply_quadrant_shifts_with_feather(
            out[j], quadrants_roi, shifts_q, margin=margin,
            shift_func=shift_func, wrap=wrap, fill=fill
        )
    return out

def plot_alignment_comparison_with_diff(
    imgs_before_cam1, imgs_before_cam2,
    imgs_after_cam1, imgs_after_cam2,
    srow, scol,
    zoom_center=None,
    zoom_size=100,
    title_prefix="Camera Alignment",
    cmap_img="gray",
    cmap_diff="seismic",
    grid_color="cyan",
    grid_step=24
):
    """
    Muestra, para cada modulación (0..3):
      - CAM1 Before (zoom)
      - CAM2 Before (zoom)
      - CAM2 - CAM1 Before (zoom)
      - CAM1 After (zoom)
      - CAM2 After (zoom)
      - CAM2 - CAM1 After (zoom)
    """
    nmods = imgs_before_cam1.shape[0]
    X, Y = imgs_before_cam1.shape[1:]

    # Centro de zoom
    cx, cy = (X // 2, Y // 2) if zoom_center is None else zoom_center

    # Ventana de recorte (clamp por seguridad)
    half = zoom_size // 2
    x0, x1 = max(0, cx - half), min(X, cx + half)
    y0, y1 = max(0, cy - half), min(Y, cy + half)

    for j in range(nmods):
        # Recortes antes
        b1 = imgs_before_cam1[j, x0:x1, y0:y1]
        b2 = imgs_before_cam2[j, x0:x1, y0:y1]
        bdiff = b2 - b1

        # Recortes después
        a1 = imgs_after_cam1[j, x0:x1, y0:y1]
        a2 = imgs_after_cam2[j, x0:x1, y0:y1]
        adiff = a2 - a1

        # Límites de visualización
        vmin_b1, vmax_b1 = _percentile_limits(b1)
        vmin_b2, vmax_b2 = _percentile_limits(b2)
        vmin_a1, vmax_a1 = _percentile_limits(a1)
        vmin_a2, vmax_a2 = _percentile_limits(a2)

        maxabs_b = np.percentile(np.abs(bdiff), 99) if bdiff.size else 1.0
        maxabs_a = np.percentile(np.abs(adiff), 99) if adiff.size else 1.0
        vmax_db = max(maxabs_b, 1e-8)
        vmax_da = max(maxabs_a, 1e-8)

        fig, axes = plt.subplots(2, 3, figsize=(13, 8))
        fig.suptitle(f"{title_prefix} – Mod {j} | shift (row, col) = ({srow:.3f}, {scol:.3f})")

        # BEFORE: CAM1
        im00 = axes[0, 0].imshow(b1, cmap=cmap_img, vmin=vmin_b1, vmax=vmax_b1)
        axes[0, 0].set_title("CAM 1 – Before")
        _apply_grid(axes[0, 0], step=grid_step, grid_color=grid_color)
        fig.colorbar(im00, ax=axes[0, 0], fraction=0.046, pad=0.04)

        # BEFORE: CAM2
        im01 = axes[0, 1].imshow(b2, cmap=cmap_img, vmin=vmin_b2, vmax=vmax_b2)
        axes[0, 1].set_title("CAM 2 – Before")
        _apply_grid(axes[0, 1], step=grid_step, grid_color=grid_color)
        fig.colorbar(im01, ax=axes[0, 1], fraction=0.046, pad=0.04)

        # BEFORE: DIFF
        im02 = axes[0, 2].imshow(bdiff, cmap=cmap_diff, vmin=-vmax_db, vmax=vmax_db)
        axes[0, 2].set_title("DIFF Before (CAM2 - CAM1)")
        _apply_grid(axes[0, 2], step=grid_step, grid_color=grid_color)
        fig.colorbar(im02, ax=axes[0, 2], fraction=0.046, pad=0.04)

        # AFTER: CAM1
        im10 = axes[1, 0].imshow(a1, cmap=cmap_img, vmin=vmin_a1, vmax=vmax_a1)
        axes[1, 0].set_title("CAM 1 – After")
        _apply_grid(axes[1, 0], step=grid_step, grid_color=grid_color)
        fig.colorbar(im10, ax=axes[1, 0], fraction=0.046, pad=0.04)

        # AFTER: CAM2
        im11 = axes[1, 1].imshow(a2, cmap=cmap_img, vmin=vmin_a2, vmax=vmax_a2)
        axes[1, 1].set_title("CAM 2 – After")
        _apply_grid(axes[1, 1], step=grid_step, grid_color=grid_color)
        fig.colorbar(im11, ax=axes[1, 1], fraction=0.046, pad=0.04)

        # AFTER: DIFF
        im12 = axes[1, 2].imshow(adiff, cmap=cmap_diff, vmin=-vmax_da, vmax=vmax_da)
        axes[1, 2].set_title("DIFF After (CAM2 - CAM1)")
        _apply_grid(axes[1, 2], step=grid_step, grid_color=grid_color)
        fig.colorbar(im12, ax=axes[1, 2], fraction=0.046, pad=0.04)

        # Métricas
        diff_before_rms = float(np.sqrt(np.mean(np.square(bdiff)))) if bdiff.size else float('nan')
        diff_after_rms  = float(np.sqrt(np.mean(np.square(adiff)))) if adiff.size else float('nan')
        diff_before_mean = float(np.mean(bdiff)) if bdiff.size else float('nan')
        diff_after_mean  = float(np.mean(adiff)) if adiff.size else float('nan')

        print(f"[Mod {j}] shift(row,col)=({srow:.3f}, {scol:.3f}) "
              f"| diff BEFORE: mean={diff_before_mean:.4g}, RMS={diff_before_rms:.4g} "
              f"| diff AFTER: mean={diff_after_mean:.4g}, RMS={diff_after_rms:.4g}")

        plt.tight_layout()
        plt.show()

def plot_two_intensity_alignment_with_diff(
    I_cam1, I_cam2_before, I_cam2_after,
    srow=0.0, scol=0.0,
    zoom_center=None, zoom_size=100,
    title_prefix="Camera alignment (demodulated I)",
    cmap_img="gray", cmap_diff="seismic",
    grid_color="cyan"
):
    """
    Compara I_cam2 antes vs después frente a I_cam1 (referencia), con zoom y diferencias.

    Parámetros
    ----------
    I_cam1 : np.ndarray (X, Y)
        Intensidad demodulada de cámara 1 (referencia).
    I_cam2_before : np.ndarray (X, Y)
        Intensidad demodulada de cámara 2 ANTES de alinear.
    I_cam2_after : np.ndarray (X, Y)
        Intensidad demodulada de cámara 2 DESPUÉS de alinear.
    srow, scol : float
        Shift aplicado (fila, columna) a cam2.
    zoom_center : (int,int) o None
        Centro (x,y) del recorte a mostrar. Si None, usa el centro de la imagen.
    zoom_size : int
        Tamaño del recorte cuadrado (px).
    title_prefix : str
        Título principal.
    cmap_img, cmap_diff : str
        Colormaps para intensidades y diferencias.
    grid_color : str
        Color de la rejilla.
    """
    assert I_cam1.ndim == 2 and I_cam2_before.ndim == 2 and I_cam2_after.ndim == 2, \
        "Las imágenes deben ser 2D."

    X, Y = I_cam1.shape

    # Definir ventana de zoom
    if zoom_center is None:
        cx, cy = X // 2, Y // 2
    else:
        cx, cy = zoom_center

    half = zoom_size // 2
    x0, x1 = max(0, cx - half), min(X, cx + half)
    y0, y1 = max(0, cy - half), min(Y, cy + half)

    # Recortes
    ref = I_cam1[x0:x1, y0:y1]
    b2  = I_cam2_before[x0:x1, y0:y1]
    a2  = I_cam2_after[x0:x1, y0:y1]

    # Diferencias
    diff_before = b2 - ref
    diff_after  = a2 - ref

    # Escalado robusto
    vmin_ref, vmax_ref = _percentile_limits(ref)
    vmin_b2,  vmax_b2  = _percentile_limits(b2)
    vmin_a2,  vmax_a2  = _percentile_limits(a2)

    vmax_db = np.percentile(np.abs(diff_before), 99)
    vmax_da = np.percentile(np.abs(diff_after),  99)
    vmax_d  = max(vmax_db, vmax_da, 1e-8)

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    fig.suptitle(f"{title_prefix}  |  shift = ({srow:.3f}, {scol:.3f})")

    # Row 0: BEFORE
    im00 = axes[0, 0].imshow(ref, cmap=cmap_img, vmin=vmin_ref, vmax=vmax_ref)
    axes[0, 0].set_title("CAM1 (I) – REF")
    axes[0, 0].grid(True, color=grid_color, linewidth=0.5, alpha=0.7)
    fig.colorbar(im00, ax=axes[0, 0], fraction=0.046, pad=0.04)

    im01 = axes[0, 1].imshow(b2, cmap=cmap_img, vmin=vmin_b2, vmax=vmax_b2)
    axes[0, 1].set_title("CAM2 (I) – BEFORE")
    axes[0, 1].grid(True, color=grid_color, linewidth=0.5, alpha=0.7)
    fig.colorbar(im01, ax=axes[0, 1], fraction=0.046, pad=0.04)

    im02 = axes[0, 2].imshow(diff_before, cmap=cmap_diff, vmin=-vmax_d, vmax=vmax_d)
    axes[0, 2].set_title("DIFF BEFORE (CAM2 - CAM1)")
    axes[0, 2].grid(True, color=grid_color, linewidth=0.5, alpha=0.7)
    fig.colorbar(im02, ax=axes[0, 2], fraction=0.046, pad=0.04)

    # Row 1: AFTER
    im10 = axes[1, 0].imshow(ref, cmap=cmap_img, vmin=vmin_ref, vmax=vmax_ref)
    axes[1, 0].set_title("CAM1 (I) – REF")
    axes[1, 0].grid(True, color=grid_color, linewidth=0.5, alpha=0.7)
    fig.colorbar(im10, ax=axes[1, 0], fraction=0.046, pad=0.04)

    im11 = axes[1, 1].imshow(a2, cmap=cmap_img, vmin=vmin_a2, vmax=vmax_a2)
    axes[1, 1].set_title("CAM2 (I) – AFTER")
    axes[1, 1].grid(True, color=grid_color, linewidth=0.5, alpha=0.7)
    fig.colorbar(im11, ax=axes[1, 1], fraction=0.046, pad=0.04)

    im12 = axes[1, 2].imshow(diff_after, cmap=cmap_diff, vmin=-vmax_d, vmax=vmax_d)
    axes[1, 2].set_title("DIFF AFTER (CAM2 - CAM1)")
    axes[1, 2].grid(True, color=grid_color, linewidth=0.5, alpha=0.7)
    fig.colorbar(im12, ax=axes[1, 2], fraction=0.046, pad=0.04)

    for ax in axes.ravel():
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
    plt.show()

    # Métricas útiles
    def _rms(x): return float(np.sqrt(np.mean(np.square(x))))
    print(f"[Camera I alignment] shift(row,col)=({srow:.3f}, {scol:.3f}) | "
          f"diff BEFORE: mean={diff_before.mean():.4g}, RMS={_rms(diff_before):.4g} | "
          f"diff AFTER:  mean={diff_after.mean():.4g},  RMS={_rms(diff_after):.4g}")
    
def compute_pseudo_images_leastsq(imgs_cam1, imgs_cam2, M1, M2, roi=None, eps=1e-12):
    """
    Construye pseudo‑imágenes 'tipo I' por modulación j, combinando CAM1 y CAM2
    con pesos α,β que:
       (i) preservan I  (α+β = 1)
      (ii) minimizan la fuga de QUV  (|| α m1 + β m2 ||)

    Además aplica una normalización fotométrica por estado (opcional) para
    equiparar las escalas de ambas cámaras antes de combinar.

    Parámetros
    ----------
    imgs_cam1, imgs_cam2 : (4, X, Y)
        Imágenes moduladas, ya alineadas entre cámaras.
    M1, M2 : (4,4)
        Matrices de modulación de cam1 y cam2. Cada fila j: [1, a, b, c].
    roi : tuple o None
        (x0, x1, y0, y1) para estimar la normalización fotométrica por estado.
        Si None, usa toda la imagen con percentiles robustos.
    eps : float
        Tolerancia numérica para evitar divisiones por cero.

    Devuelve
    --------
    pseudo : (4, X, Y)
        Pseudo‑imágenes por modulación (aprox. 'I').
    alphas, betas : (4,)
        Pesos por estado j.
    gammas : (4,)
        Factores de normalización aplicados a cam2 por estado j.
    residual_norms : (4,)
        Norma de la fuga residual || α m1 + β m2 || por estado j.
    """
    nmods, X, Y = imgs_cam1.shape
    pseudo = np.zeros_like(imgs_cam1, dtype=np.float64)
    alphas = np.zeros(nmods)
    betas  = np.zeros(nmods)
    gammas = np.ones(nmods)
    residual_norms = np.zeros(nmods)

    for j in range(nmods):
        I1 = imgs_cam1[j].astype(np.float64)
        I2 = imgs_cam2[j].astype(np.float64)

        # --- Normalización fotométrica por estado j (cam2 -> cam1) ---
        if roi is not None:
            x0, x1, y0, y1 = roi
            ref1 = I1[x0:x1, y0:y1]
            ref2 = I2[x0:x1, y0:y1]
        else:
            # Percentiles robustos en toda la imagen
            ref1 = I1
            ref2 = I2

        # Evitar outliers (usar median de ratios donde I1 no es casi cero)
        mask = np.abs(ref1) > np.percentile(np.abs(ref1), 5) + eps
        if np.any(mask):
            ratio = ref2[mask] / (ref1[mask] + eps)
            gamma = np.median(ratio)
        else:
            gamma = 1.0

        I2n = I2 / (gamma + eps)
        gammas[j] = gamma

        # --- Pesos óptimos por mínimos cuadrados con α+β=1 ---
        m1 = M1[j, 1:].astype(np.float64)  # (a,b,c)
        m2 = M2[j, 1:].astype(np.float64)  # (a',b',c')

        diff = m1 - m2
        denom = np.dot(diff, diff)
        if denom < eps:
            # Caso degenerado: m1 ≈ m2 → cualquier α con α+β=1 sirve; usa α=β=0.5
            alpha = 0.5
        else:
            alpha = (np.dot(m2, m2) - np.dot(m1, m2)) / denom
        beta = 1.0 - alpha

        alphas[j] = alpha
        betas[j]  = beta

        # Fuga residual teórica (diagnóstico)
        residual = alpha * m1 + beta * m2
        residual_norms[j] = np.linalg.norm(residual)

        # --- Pseudo‑imagen ---
        pseudo[j] = alpha * I1 + beta * I2n

        # (Opcional) Re‑escala para visualizar comparables
        # pseudo[j] *= 1.0  # mantener unidad de I por construcción (α+β=1)

    return pseudo, alphas, betas, gammas, residual_norms

def plot_pseudo_mod_alignment(
    pseudo_before, pseudo_after,
    srow_mod, scol_mod,
    zoom_center=None,
    zoom_size=100,
    title_prefix="Pseudo‑I modulation alignment",
    cmap_img="gray",
    cmap_diff="seismic",
    grid_color="yellow"
):
    """
    Compara pseudo‑imágenes de cada modulación j (1,2,3) frente a la mod 0,
    mostrando imágenes y diferencias antes y después del alineamiento.

    pseudo_before : (4, X, Y)
    pseudo_after  : (4, X, Y)
    """
    X, Y = pseudo_before.shape[1:]

    # --- ZOOM window ---
    if zoom_center is None:
        cx, cy = X//2, Y//2
    else:
        cx, cy = zoom_center
    half = zoom_size//2
    x0, x1 = max(0, cx-half), min(X, cx+half)
    y0, y1 = max(0, cy-half), min(Y, cy+half)

    # --- Referencia mod 0 ---
    ref_b = pseudo_before[0, x0:x1, y0:y1]
    ref_a = pseudo_after[0,  x0:x1, y0:y1]

    for j in [1, 2, 3]:
        img_b = pseudo_before[j, x0:x1, y0:y1]
        img_a = pseudo_after[j,  x0:x1, y0:y1]

        diff_b = img_b - ref_b
        diff_a = img_a - ref_a

        # Límites
        vmin_r, vmax_r = np.percentile(ref_b, 1), np.percentile(ref_b, 99)
        vmin_j, vmax_j = np.percentile(img_b, 1), np.percentile(img_b, 99)

        vmax_diff_b = np.percentile(np.abs(diff_b), 99)
        vmax_diff_a = np.percentile(np.abs(diff_a), 99)

        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        fig.suptitle(f"{title_prefix} – Compare mod 0 vs mod {j} "
                     f" | shift({j})=({srow_mod[j]:.3f}, {scol_mod[j]:.3f})")

        # BEFORE: ref mod 0
        im00 = axes[0, 0].imshow(ref_b, cmap=cmap_img, vmin=vmin_r, vmax=vmax_r)
        axes[0, 0].set_title("Mod 0 BEFORE")
        axes[0, 0].grid(True, color=grid_color, linewidth=0.5)
        fig.colorbar(im00, ax=axes[0, 0])

        # BEFORE: mod j
        im01 = axes[0, 1].imshow(img_b, cmap=cmap_img, vmin=vmin_j, vmax=vmax_j)
        axes[0, 1].set_title(f"Mod {j} BEFORE")
        axes[0, 1].grid(True, color=grid_color, linewidth=0.5)
        fig.colorbar(im01, ax=axes[0, 1])

        # BEFORE: difference
        im02 = axes[0, 2].imshow(diff_b, cmap=cmap_diff, vmin=-vmax_diff_b, vmax=vmax_diff_b)
        axes[0, 2].set_title(f"Diff BEFORE (mod {j} - mod 0)")
        axes[0, 2].grid(True, color=grid_color, linewidth=0.5)
        fig.colorbar(im02, ax=axes[0, 2])

        # AFTER: ref mod 0
        im10 = axes[1, 0].imshow(ref_a, cmap=cmap_img, vmin=vmin_r, vmax=vmax_r)
        axes[1, 0].set_title("Mod 0 AFTER")
        axes[1, 0].grid(True, color=grid_color, linewidth=0.5)
        fig.colorbar(im10, ax=axes[1, 0])

        # AFTER: mod j
        im11 = axes[1, 1].imshow(img_a, cmap=cmap_img, vmin=vmin_j, vmax=vmax_j)
        axes[1, 1].set_title(f"Mod {j} AFTER")
        axes[1, 1].grid(True, color=grid_color, linewidth=0.5)
        fig.colorbar(im11, ax=axes[1, 1])

        # AFTER: difference
        im12 = axes[1, 2].imshow(diff_a, cmap=cmap_diff, vmin=-vmax_diff_a, vmax=vmax_diff_a)
        axes[1, 2].set_title(f"Diff AFTER (mod {j} - mod 0)")
        axes[1, 2].grid(True, color=grid_color, linewidth=0.5)
        fig.colorbar(im12, ax=axes[1, 2])

        plt.tight_layout()
        plt.show()

        # Métricas
        print(f"[Pseudo] mod 0 vs mod {j}:")
        print(f"    BEFORE diff RMS = {np.sqrt(np.mean(diff_b**2)):.4g}")
        print(f"    AFTER  diff RMS = {np.sqrt(np.mean(diff_a**2)):.4g}")

def dftreg(F,G,kappa):
    """
    Calculates the shift between a couple of images 'f' and 'g' with subpixel
    accuracy following the second method presented in
    Sicairos 2008, Efficient subpixel image registration algorithm.
    Input:
        F,G: ffts of images 'f' and 'g' without applying any fftshift
        kappa: inverse of subpixel precision (kappa=20 -> 0.05 pixel precision)
    Output:

    """
    nr,nc=np.shape(F)
    Nr = ifftshift(np.arange(-np.fix(nr/2),np.ceil(nr/2)))
    Nc = ifftshift(np.arange(-np.fix(nc/2),np.ceil(nc/2)))
    Nout = 2 * max(nr, nc)
    CC=ifft2(FTpad(F*np.conj(G), Nout))
    CCabs=np.abs(CC)
    ind = np.unravel_index(np.argmax(CCabs, axis=None), CCabs.shape)

    CCmax=CC[ind]*nr*nc
    Nr2 = ifftshift(np.arange(-np.fix(nr),np.ceil(nr)))
    Nc2 = ifftshift(np.arange(-np.fix(nc),np.ceil(nc)))

    row_shift=Nr2[ind[0]]/2
    col_shift=Nc2[ind[1]]/2

    #Initial shift estimate in upsampled grid
    row_shift=round(row_shift*kappa)/kappa
    col_shift=round(col_shift*kappa)/kappa
    dftshift=np.fix(np.ceil(kappa*1.5)/2)

    #DFT by matrix multiplication
    CC=np.conj(dftups(G*np.conj(F),np.ceil(kappa*1.5),kappa,\
    dftshift-row_shift*kappa,dftshift-col_shift*kappa))
    CCabs=np.abs(CC)
    ind = np.unravel_index(np.argmax(CCabs, axis=None), CCabs.shape)
    CCmax=CC[ind]
    rloc,cloc=ind-dftshift
    row_shift=row_shift+rloc/kappa
    col_shift=col_shift+cloc/kappa
    rg00=np.sum(np.abs(F)**2)
    rf00=np.sum(np.abs(G)**2)
    error=np.sqrt(1-np.abs(CCmax)**2/(rg00*rf00))
    Nc,Nr=np.meshgrid(Nc,Nr)
    Gshift=G*np.exp(1j*2*np.pi*(-row_shift*Nr/nr-col_shift*Nc/nc))
    return error,row_shift,col_shift,Gshift

def dftups(M,n_out,kappa,roff,coff):
    """
    Upsampled cross-correlation obtained by matrix multiplication
    Inputs:
        M: input image for calculation of the DFT
        n_out: number of pixels in the output upsampled DFT
        kappa: inverse of subpixel precision (kappa=20 -> 0.005 pixel precision)
        roff, coff: row and column offsets to shift the output array to a
            region of interest
    """
    nr,nc=M.shape
    kernc=np.exp((-1j*2*np.pi/(nc*kappa))*np.outer(\
    ifftshift(np.arange(0,nc).T-np.floor(nc/2)),np.arange(0,n_out)-coff))

    kernr=np.exp((-1j*2*np.pi/(nr*kappa))*np.outer(\
    np.arange(0,n_out)-roff,ifftshift(np.arange(0,nr).T-np.floor(nr/2))))
    return kernr @ M @ kernc

def FTpad(IM,Nout):
    """
    Carries out zero-padding to upsample an image IM in Fourier domain
    Input:
        IM: Numpy array in Fourier domain
        outsize: size of the new array

    """
    Nin=IM.shape[0]
    pd=int((Nout-Nin)/2)
    IM=fftshift(IM)
    IMout=np.pad(IM,((pd,pd),(pd,pd)),'constant')
    IMout=ifftshift(IMout)*Nout*Nout/(Nin*Nin)
    return IMout

# ------------------------------  MAIN FUNCTS  --------------------------------- # 

def realign_subpixel(ima, accu=0.01, verbose = True, return_shift = False):
    """
    This function aligns a series of images with subpixel images using the Sicairos
    method.
    Input:
     ima: 3D array of the type (Nima, Nx, Ny). First dimension corresponds to the
        index of the image through the series
     accu: accuracy of the alignment in pixel units
    Output: returns the aligned 3D array
    """
    kappa = 1 / accu #Kappa factor defined in Sicairos method (1/fraction of pixel)
    Gshift = fft2(ima[0, :, :]) # FFT of the first image of the series
    if verbose:
        print('Re-aligning images ...')  
    ima_aligned = np.zeros(np.shape(ima))
    
    row_shifts = []
    col_shifts = []
    for j in range(ima.shape[0]):
        
        F0=np.copy(Gshift)

        F_comp = fft2(ima[j])
        error, row_shift, col_shift, Gshift2 = dftreg(F0, F_comp, kappa)
        row_shifts.append(row_shift)
        col_shifts.append(col_shift)
        if verbose:
            print(f"Shift of image: {j} -> row : {round(row_shift, 4)} col : {round(col_shift, 4)}")
        if j != 0:
            ima_aligned[j] = np.real(ifft2(Gshift2))
        else:
            ima_aligned[j] = ima[0]
    
    if return_shift:
        return ima_aligned, row_shifts, col_shifts, error
    else:
        return ima_aligned, error

def find_fieldstop(cam1 = None, verbose = False, plot_flag = False, margin = 10):
    """
    Module to find the fieldstop of images. 

    Inputs:
        - cam1 (np.array): A single image to find the fieldstop.
        - verbose (Boolean, default : False): Print info on terminal. 
        - plot_flag (Boolean, default : False): Plot the fieldstop calculation. 
        - margin (int, default : 10) : Number of pixels of margin from the detected field-stop       
    Outputs:
        - Fieldstop (list). 
    """

    tic = time.time()

    if verbose:
        print("Finding fieldstop field stop...")
    
    if plot_flag:
        fig, axs  = plt.subplots(1, 2,figsize = (10, 5))
        axs[0].imshow(cam1, origin = 'lower', cmap = 'gray')
    
    size = np.shape(cam1)[0]

    # Position to find cuts
    lines = np.linspace(0, size, 7)
    lines = [int(x) for x in lines[1:-1]]
    
    # Looking for cuts
    hcuts_left_c1 = []
    hcuts_right_c1 = []
    vcuts_top_c1 = []
    vcuts_bottom_c1 = []
    for l in lines:
        # Camera 1
        hcut1 = np.argmax(np.gradient(cam1[l, :]))
        hcut2 = np.argmin(np.gradient(cam1[l, :]))
        vcut1 = np.argmax(np.gradient(cam1[:, l]))
        vcut2 = np.argmin(np.gradient(cam1[:, l]))
        hcuts_right_c1.append(hcut1)
        hcuts_left_c1.append(hcut2)                                    
        vcuts_top_c1.append(vcut1)
        vcuts_bottom_c1.append(vcut2)

        if plot_flag:
            axs[0].plot([l, l], [0, size], color = 'crimson' , lw = 1)
            axs[0].plot([0, size], [l, l], color = 'crimson' , lw = 1)
            axs[0].scatter(l, vcut1, marker = 'x', c = 'dodgerblue')
            axs[0].scatter(l, vcut2, marker = 'x', c = 'darkorange')
            axs[0].scatter(hcut1, l, marker = 'x', c = 'dodgerblue')
            axs[0].scatter(hcut2, l, marker = 'x', c = 'darkorange')

    # Selecting the innermost points (in case border is tilted)

    vcut_right_c1 = np.min(hcuts_left_c1) - margin 
    vcut_left_c1 = np.max(hcuts_right_c1) + margin
    hcut_top_c1 = np.min(vcuts_bottom_c1) - margin
    hcut_bottom_c1 = np.max(vcuts_top_c1) + margin

    cam1_fieldstop = np.array([[hcut_bottom_c1, hcut_top_c1], [vcut_left_c1, vcut_right_c1]])

    if plot_flag:
        axs[0].plot([vcut_right_c1, vcut_right_c1], [0, size], c = 'deeppink')
        axs[0].plot([vcut_left_c1, vcut_left_c1], [0, size], c = 'deeppink')
        axs[0].plot([0, size], [hcut_top_c1, hcut_top_c1], c = 'deeppink')
        axs[0].plot([0, size], [hcut_bottom_c1, hcut_bottom_c1], c = 'deeppink')
        axs[1].imshow(cam1[hcut_bottom_c1:hcut_top_c1, vcut_left_c1:vcut_right_c1], origin = 'lower', cmap = 'gray')
        axs[0].set_xlim(0, size)
        axs[0].set_ylim(0, size)
        axs[0].set_ylabel("Cam 1")
        plt.tight_layout()
        plt.show()

    print(f"Field stop computation finished in {round(time.time() - tic, 3)}s.")

    return cam1_fieldstop
  
def shift_subp(im: np.ndarray, shift=None, wrap=True, fill=0):
    '''define shift operator (subpixel)
        Input is y and x shifts (defined negative towards (0,0)
        new center = center + (x,y)
        Note that image is defined as [sy,sx] so shifts = [sy (rows),sx (columns)]
    '''
    import math
    nr, nc = im.shape
    Nr = ifftshift(np.arange(-np.fix(nr / 2), np.ceil(nr / 2)))
    Nc = ifftshift(np.arange(-np.fix(nc / 2), np.ceil(nc / 2)))
    Nc, Nr = np.meshgrid(Nc, Nr)
    G = fft2(im)
    Gshift = G * np.exp(1j * 2 * np.pi * (-shift[0] * Nr / nr - shift[1] * Nc / nc))
    im_shift = np.real(ifft2(Gshift))

    if wrap is False:
        dy, dx = shift
        if dx > 0:
            im_shift[:, 0:math.ceil(dx)] = fill
        elif dx < 0:
            im_shift[:, math.floor(dx):] = fill
        if dy > 0:
            im_shift[0:math.ceil(dy), :] = fill
        elif dy < 0:
            im_shift[math.floor(dy):, :] = fill

    return im_shift

def align_obsmode(data, acc = 0.01, verbose = False, filter = filter, 
                  onelambda = False, returnshifts = True,roi = [0,-1,0,-1], quadrants = 0,
                  align_sequence = 0, debug = False):
    """
    Function to filter, rotate camera 2 and align an obs mode. 

    Inputs:
        - data (np.array) : Array contaning the obs mode. (Ncams x Nlambda x Nmods x Nx x Ny)
        - acc (float,default : 0.01) : Accuracy for the alignemnt routine.
        - theta (float, default : 0.0655): Angle of rotation
        - verbose (Boolean, ddefault : False) : Print info on terminal.
        - filterflag (Boolean, default : True) : Set to False to skip Fourier filtration 
        - onelambda (Boolean, default : False): Set to true if only one lambda is used (array of shape Ncam x Nmod x Nx x Ny) 
        - zkes (np.array, default : np.zeros(21)): Zernike's array to use for the filtration. 
    Outputs:
        - Filtered, Rotated  and aligned (np.array) : Same array as data filtrated and with cam2 rotated  
        - shifts (list) : shifts performed to each camera, modulation and wavelength

    Notas: Sobre la opción modo = 1. Solo activa si se activa. Jaja

    Procesamiento para alineación en TuMag con dos cámaras y modulaciones ortogonales.
    =================================================================================

    Esta opción implementa el pipeline de alineación suponiendo (en primer orden)
    que las dos cámaras de TuMag registran *el mismo jitter temporal*. Sin embargo,
    las cámaras miden modulaciones *ortogonales*, de modo que las imágenes moduladas
    de una cámara y de la otra no son directamente comparables entre sí. Tampoco lo
    son dentro de una misma cámara, ya que cada imagen corresponde a una modulación
    distinta.

    Descripción del problema
    ------------------------
    Cada cámara registra cuatro estados de modulación distintos:

        I_j^(cam1) = I + M^(1)_{1j} Q + M^(1)_{2j} U + M^(1)_{3j} V
        I_j^(cam2) = I + M^(2)_{1j} Q + M^(2)_{2j} U + M^(2)_{3j} V

    donde cada conjunto M^(1) y M^(2) corresponde a modulaciones ortogonales, y
    el índice j identifica el estado de modulación.

    Aunque ambas cámaras sufren exactamente el mismo jitter temporal, sus imágenes
    no pueden alinearse estado a estado, ya que cada una es una combinación lineal
    distinta de los parámetros de Stokes.

    Objetivo
    --------
    Obtener los parámetros de Stokes (I, Q, U, V) correctamente alineados entre
    cámaras, eliminando los efectos del jitter y garantizando coherencia geométrica
    tanto a nivel de las modulaciones como en los Stokes finales demodulados.

    Método
    ------
    El procedimiento (idealizado) implementado es el siguiente:

    1. **Demodulación independiente por cámara.**
    Se calculan (I, Q, U, V) de cada cámara por separado sin intentar alinear
    aún las modulaciones crudas.

    2. **Alineación utilizando únicamente Stokes I.**
    Puesto que Stokes I debe ser idéntico en ambas cámaras salvo el desplazamiento
    relativo entre ellas, se calcula la transformación geométrica (desplazamientos)
    que satisface:
        I_cam1 → I_cam2
        I_cam2 → I_cam1
    Para evitar sesgos, se promedian ambas transformaciones, obteniendo una
    transformación común, apropiada para aplicarse a las modulaciones crudas.

    3. **Aplicación de la transformación obtenida a las modulaciones crudas.**
    La transformación geométrica derivada de Stokes I se aplica ahora a las
    imágenes moduladas originales, antes de la demodulación final. Es suficiente
    aplicarla a una sola cámara; por comodidad se aplica a la cámara 2.

    4. **Alineación fina de las modulaciones.**
    Debido a que cada estado de modulación es una combinación lineal distinta de
    (I, Q, U, V), para poder alinear los estados correctamente se generan
    combinaciones que atenúan Q, U y V, obteniendo pseudo‑intensidades comparables
    entre pseudo‑modulaciones. Estas señales permiten determinar con mayor
    precisión el jitter relativo entre modulaciones, bajo la suposición inicial
    de que dicho jitter es idéntico en ambas cámaras.

    5. **Corrección final de modulaciones.**
    Una vez obtenidos los desplazamientos entre modulaciones, se corrigen todos
    los estados en ambas cámaras utilizando el mismo jitter. Es necesario tener
    en cuenta la orientación relativa de las cámaras para aplicar correctamente
    estas correcciones.

    Consideraciones adicionales
    ---------------------------
    - En principio, la alineación debería ser casi constante con la longitud de
    onda, pero efectos ópticos (rotaciones, distorsiones, filtros, etc.)
    introducen pequeñas variaciones con longitud de, por lo que no se fija una
    transformación global entre cámaras. Se evalúa por cada landa.

    """

    def _generar_cuadrantes_nx_n(H, W, n):
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

    tic = time.time() # Get the time to measure execution time.

    if onelambda:
        data = data[:, np.newaxis] # To allow for only one lamdba.

    shape = np.shape(data)
    nlambda = shape[1]
    nmods = shape[2]

    aligned =  np.copy(data)

    if quadrants != 0:
        shifts = np.zeros((nlambda, abs(quadrants*quadrants), 2, 2, nmods), dtype=float)
        # ojo, el ROI tiene que coincidir con la imagen así que lo fuerzo
        roi[0] = 0
        roi[1] = data.shape[-2]
        roi[2] = 0
        roi[3] = data.shape[-1]

        if quadrants > 0:
            H, W = data.shape[-1], data.shape[-1]   # dimensiones de la imagen
            quadrants_roi = _generar_cuadrantes_nx_n(H, W, quadrants)
            # for idx, (y1, y2, x1, x2) in enumerate(quadrants_roi):
            #     print(f"ROI {idx}: y={y1}:{y2}, x={x1}:{x2}")
    else:
        shifts = np.zeros((nlambda, 2, 2, nmods))

    if verbose:
        plt.imshow(data[0,0,0,roi[0]:roi[1],roi[2]:roi[3]])
        plt.show()
    err = []

    for lambd in range(nlambda):
        logging.info(f"Aligning wavelength: {lambd + 1}/{nlambda}")

        if quadrants == 0:
            logging.info(f"quadrants 0: {quadrants}")
            if align_sequence == 0:
                if verbose:
                    logging.info(f"Shifts for cam 1 - modulation alignment")

                _, srow, scol, error = realign_subpixel(data[0, lambd,:,roi[0]:roi[1],roi[2]:roi[3]], verbose = verbose, accu = acc, return_shift=True)
                err.append(error)

                shifts[lambd, 0, 0] = srow
                shifts[lambd, 0, 1] = scol

                for nm in range(nmods-1):
                    aligned[0, lambd,nm + 1] = shift_subp(data[0, lambd,nm+1], shift=[srow[nm + 1], scol[nm + 1]], wrap=True, fill=0)

                if verbose:
                    logging.info("Shifts of camera 2 alignment")
                for mod in range(nmods):

                    _, srow, scol,error = realign_subpixel(np.array([aligned[0,lambd,mod,roi[0]:roi[1],roi[2]:roi[3]], data[1, lambd, mod,roi[0]:roi[1],roi[2]:roi[3]]]), verbose = verbose, accu = acc, return_shift=True)
                    err.append(error)

                    shifts[lambd, 1, 0, mod] = srow[1]
                    shifts[lambd, 1, 1, mod] = scol[1]

                    aligned[1, lambd,mod] = shift_subp(data[1, lambd,mod], shift=[srow[1], scol[1]], wrap=True, fill=0)

            elif align_sequence == 1:

                # =============================
                # (1) DEMODULACIÓN DE LAS CAMARAS
                # =============================
                _, demod = demodulate(data[:, lambd], 
                                    filt=filter, 
                                    # dmod_matrices='demod_matrices_david',
                                    onelambda=True, 
                                    BothCams=True)

                # demod tiene shape (2, 4, x, y)
                I_cam1 = demod[0, 0]
                I_cam2 = demod[1, 0]

                scale, gamma = balance(I_cam1, I_cam2, roi=roi)
                I_cam2 = I_cam2 * scale
                logging.info(f"Balance (2D): scale={scale:.6f}, gamma={gamma:.6f}")

                # =============================
                # (2) ALINEAR CAM1 <--> CAM2 CON STOKES I
                # =============================
                logging.info("Global alignment between cameras using Stokes I")

                # TODO, Falta hacerlo al revés y promediar
                _, srow_cam, scol_cam, err_I = realign_subpixel(
                    np.array([I_cam1[roi[0]:roi[1], roi[2]:roi[3]],
                            I_cam2[roi[0]:roi[1], roi[2]:roi[3]]]),
                    verbose=verbose, accu=acc, return_shift=True
                )
                # Print numeric shifts
                logging.info(f"Global camera shift → row={srow_cam[1]:.3f}, col={scol_cam[1]:.3f}")

                # _, srow_cam, scol_cam, err_I = realign_subpixel(
                #     np.array([I_cam2[roi[0]:roi[1], roi[2]:roi[3]],
                #             I_cam1[roi[0]:roi[1], roi[2]:roi[3]]]),
                #     verbose=verbose, accu=acc, return_shift=True
                # )
                # # Print numeric shifts
                # print(f"Global camera shift → row={-srow_cam[1]:.3f}, col={-scol_cam[1]:.3f}")

                err.append(err_I)

                # Guardamos shifts globales
                shifts[lambd, 1, 0, :] = srow_cam[1]   # cam2 shift
                shifts[lambd, 1, 1, :] = scol_cam[1]

                # Aplicar shift global SOLO a cam2 
                for j in range(nmods):
                    aligned[1, lambd, j] = shift_subp(
                        data[1, lambd, j],
                        shift=[srow_cam[1], scol_cam[1]],
                        wrap=True, fill=0
                    )
                # Cam1 no se toca (referencia)
                aligned[0, lambd] = data[0, lambd]

                if debug['activate']:

                    I_cam2_before = I_cam2.copy()

                    I_cam2 = shift_subp(
                        I_cam2,
                        shift=[srow_cam[1], scol_cam[1]],
                        wrap=True,
                        fill=0
                    )

                    # Visualización
                    plot_two_intensity_alignment_with_diff(
                        I_cam1=I_cam1,
                        I_cam2_before=I_cam2_before,
                        I_cam2_after=I_cam2,
                        srow=srow_cam[1],
                        scol=scol_cam[1],
                        zoom_center=debug['zoom_center'],  # p.ej. (400, 400)
                        zoom_size=debug['zoom_size'],      # p.ej. 100
                        title_prefix=f"λ={lambd} Camera alignment before vs after"
                    )

                    # ---------- DEBUG VISUALIZATION ----------
                    # Show before/after alignment for this λ

                    imgs_before_cam1 = data[0, lambd]           # (4,x,y)
                    imgs_before_cam2 = data[1, lambd]           # (4,x,y)

                    imgs_after_cam1 = aligned[0, lambd]         # cam1 reference
                    imgs_after_cam2 = aligned[1, lambd]         # cam2 shifted


                    # Display comparison
                    plot_alignment_comparison_with_diff(
                        imgs_before_cam1, imgs_before_cam2,
                        imgs_after_cam1, imgs_after_cam2,
                        srow=srow_cam[1],
                        scol=scol_cam[1],
                        zoom_center=debug['zoom_center'],
                        zoom_size=debug['zoom_size'],
                        title_prefix=f"λ={lambd} Camera alignment (M1_cam1 M1_cam2)"
                    )
                # -------------------------------------------


                # =============================
                # (3) PSEUDO-IMÁGENES POR MODULACIÓN j
                # =============================
                imgs_cam1 = aligned[0, lambd]    # (4,x,y)
                imgs_cam2 = aligned[1, lambd]    # (4,x,y)

                M1 = mod_matrices[filter][0]
                M2 = mod_matrices[filter][1]


                roi_norm = (roi[0], roi[1], roi[2], roi[3])

                pseudo_imgs, alphas, betas, gammas, residuals = compute_pseudo_images_leastsq(
                    imgs_cam1, imgs_cam2, M1, M2, roi=roi_norm
                )

                for j in range(4):
                    print(f"[λ={lambd}] mod {j}: alpha={alphas[j]:.4f}, beta={betas[j]:.4f}, "
                        f"gamma={gammas[j]:.4f}, residual_norm={residuals[j]:.4e}")

                # =============================
                # (4) ALINEACIÓN ENTRE MODULACIONES
                # =============================
                logging.info("Aligning modulations using pseudo-images...")

                _, srow_mod, scol_mod, err_mod = realign_subpixel(
                    pseudo_imgs[:, roi[0]:roi[1], roi[2]:roi[3]],
                    verbose=verbose, accu=acc, return_shift=True
                )
                err.append(err_mod)

                shifts[lambd, 0, 0] = srow_mod
                shifts[lambd, 0, 1] = scol_mod

                for j in range(4):
                    print(f"Camera shift mod {j} → row={shifts[lambd, 0, 0][j]:.3f}, col={shifts[lambd, 0, 1][j]:.3f}")

                if debug['activate']:
                    # ===== DEBUG: Visualización del alineamiento de pseudo‑imágenes =====

                    pseudo_before = pseudo_imgs.copy()

                    # Construir pseudo_after aplicando el jitter
                    pseudo_after = np.zeros_like(pseudo_imgs)
                    for j in range(nmods):
                        pseudo_after[j] = shift_subp(
                            pseudo_imgs[j],
                            shift=[srow_mod[j], scol_mod[j]],
                            wrap=True,
                            fill=0
                        )

                    plot_pseudo_mod_alignment(
                        pseudo_before,
                        pseudo_after,
                        srow_mod, scol_mod,
                        zoom_center=debug['zoom_center'],
                        zoom_size=debug['zoom_size'],
                        title_prefix=f"λ={lambd} – pseudo‑modulation alignment"
                    )


                # =============================
                # (5) APLICAR ESTE JITTER A AMBAS CAMARAS
                # =============================
                for cam in range(2):
                    for j in range(nmods):
                        aligned[cam, lambd, j] = shift_subp(
                            aligned[cam, lambd, j],
                            shift=[srow_mod[j], scol_mod[j]],
                            wrap=True,
                            fill=0
                        )

                logging.info("Final modulation alignment applied.")

            else:
                exit

        elif quadrants > 0:
            logging.info(f"quadrants: {quadrants}")

            if align_sequence == 1:
                # =============================
                # (1) DEMODULACIÓN DE LAS CAMARAS (sin BothCams)
                # =============================
                _, demod = demodulate(
                    data[:, lambd],
                    filt=filter,
                    # dmod_matrices='demod_matrices_david',
                    onelambda=True,
                    BothCams=True
                )
                # demod: (2, 4, x, y)
                I_cam1 = demod[0, 0].astype(float)  # I de cam1
                I_cam2 = demod[1, 0].astype(float)  # I de cam2

                # =============================
                # (1b) BALANCE (2D) SOLO EN ESTA λ (opcional pero recomendado)
                # =============================
                scale, gamma = balance(I_cam1, I_cam2, roi=roi)
                I_cam2_bal = I_cam2 * scale
                print(f"[λ={lambd}] Balance (2D): scale={scale:.6f}, gamma={gamma:.6f}")

                # =============================
                # (2) ALINEAR CAM1 <--> CAM2 CON STOKES I (POR CUADRANTE)
                # =============================
                logging.info("Camera alignment (per quadrant) using Stokes I")

                # Empezamos con 'aligned' como copia de los datos modulados crudos
                aligned[0, lambd] = data[0, lambd]  # cam1 como referencia
                aligned[1, lambd] = data[1, lambd]  # luego sustituiremos cuadrantes de cam2

                for q, (y1, y2, x1, x2) in enumerate(quadrants_roi):

                    patch1 = I_cam1[y1:y2, x1:x2]
                    patch2 = I_cam2_bal[y1:y2, x1:x2]

                    # Shift cam2 -> cam1 en ESTE cuadrante
                    _, srow_cam_q, scol_cam_q, _ = realign_subpixel(
                        np.array([patch1, patch2]),
                        verbose=verbose, accu=acc, return_shift=True
                    )
                    dy_cam = float(srow_cam_q[1])
                    dx_cam = float(scol_cam_q[1])

                    # Guardar shifts por cuadrante para "cámara 2"
                    # Estructura: shifts[lambda, q, cam, axis, mod]
                    # Para cámara 2 (cam=1), rellenamos el mismo shift para los 4 mods
                    shifts[lambd, q, 1, 0, :] = dy_cam
                    shifts[lambd, q, 1, 1, :] = dx_cam

                    print(f"  Q{q}: cam2 shift (dy,dx)=({dy_cam:.3f}, {dx_cam:.3f})")

                # shifts por cuadrante para CAM2 (cam2 -> cam1) — vector de Q elementos
                shifts_cam_q = [(float(shifts[lambd, q, 1, 0, 0]), float(shifts[lambd, q, 1, 1, 0]))
                                for q in range(len(quadrants_roi))]

                # Aplicar el shift de cámara 2 SOLO en este cuadrante a los 4 estados modulados
                for j in range(nmods):
                    # aligned[1, lambd, j] = _apply_patch_shift(
                    #     aligned[1, lambd, j], y1, y2, x1, x2, dy_cam, dx_cam,
                    #     wrap=True, fill=0
                    # )
                    # TODO: SIGUE HABIENDO UN PROBLEMA AL UNIR LAS IMAGENES
                    aligned[1, lambd, j] = apply_quadrant_shifts_with_feather(
                        aligned[1, lambd, j],
                        quadrants_roi=quadrants_roi,
                        shifts_q=shifts_cam_q,    # mismos shifts de cámara para todos los mods
                        margin=10,
                        shift_func=shift_subp,
                        wrap=False, fill=0
                    )
                
                # =============================
                # (2c) DEBUG: DEMODULADAS (I) POR CUADRANTE (OPCIONAL)
                # =============================
                if debug['activate']:

                    H, W = shape[-1], shape[-2]  # pon aquí tus valores
                    dy_cam2, dx_cam2 = build_displacement_maps_simple(H, W, quadrants_roi, shifts_cam_q)
                    plot_maps(dy_cam2, dx_cam2, title="CAM2")

                    # Usaremos la imagen balanceada como "before" y aplicaremos los shifts por cuadrante sobre ella.
                    I_cam2_before = I_cam2_bal.copy()

                    # Shifts por cuadrante para CAM2 (cam2 -> cam1)
                    # Nota: hemos guardado (para cam=1) el mismo shift de cámara en los 4 mods; cogemos el de mod 0.
                    shifts_cam_q = []
                    for q, (y1, y2, x1, x2) in enumerate(quadrants_roi):
                        dy = float(shifts[lambd, q, 1, 0, 0])  # row shift, cam2, mod 0
                        dx = float(shifts[lambd, q, 1, 1, 0])  # col shift, cam2, mod 0
                        shifts_cam_q.append((dy, dx))
                        if verbose:
                            print(f"  [DEBUG] λ={lambd} Q{q}: cam2 shift (dy,dx)=({dy:+.3f}, {dx:+.3f})")

                    # Aplicar shifts por cuadrante con apodización y mezcla suave
                    I_cam2_after = apply_quadrant_shifts_with_feather(
                        img2d=I_cam2_before,
                        quadrants_roi=quadrants_roi,
                        shifts_q=shifts_cam_q,
                        margin=10,
                        shift_func=shift_subp,
                        wrap=True,   # importante para no sangrar contenido de otros cuadrantes
                        fill=0
                    )

                    # Visualización comparativa (antes vs después) frente a I_cam1
                    # Aquí no mostramos un único (srow, scol) porque trabajamos por cuadrante
                    plot_two_intensity_alignment_with_diff(
                        I_cam1=I_cam1,
                        I_cam2_before=I_cam2_before,
                        I_cam2_after=I_cam2_after,
                        srow=0.0, scol=0.0,  # per-quadrant → no aplica un único valor
                        zoom_center=debug.get('zoom_center', (400, 400)),
                        zoom_size=debug.get('zoom_size', 100),
                        title_prefix=f"λ={lambd} Camera alignment (per quadrant, feather={10}px)"
                    )

                # =============================
                # (3) PSEUDO-IMÁGENES POR MODULACIÓN j (con cámaras ya patch‑alineadas)
                # =============================
                imgs_cam1 = aligned[0, lambd]    # (4,x,y)
                imgs_cam2 = aligned[1, lambd]    # (4,x,y)

                M1 = mod_matrices[filter][0]
                M2 = mod_matrices[filter][1]

                roi_norm = (roi[0], roi[1], roi[2], roi[3])
                pseudo_imgs, alphas, betas, gammas, residuals = compute_pseudo_images_leastsq(
                    imgs_cam1, imgs_cam2, M1, M2, roi=roi_norm
                )

                for j in range(nmods):
                    print(f"[λ={lambd}] mod {j}: alpha={alphas[j]:.4f}, beta={betas[j]:.4f}, "
                        f"gamma={gammas[j]:.4f}, residual_norm={residuals[j]:.4e}")

                # =============================
                # (4) ALINEACIÓN ENTRE MODULACIONES (POR CUADRANTE)
                # =============================
                logging.info("Modulation alignment (per quadrant) using pseudo‑images")

                for q, (y1, y2, x1, x2) in enumerate(quadrants_roi):
                    patch_pseudo = pseudo_imgs[:, y1:y2, x1:x2]  # (4, h, w)
                    _, srow_mod_q, scol_mod_q, _ = realign_subpixel(
                        patch_pseudo, verbose=verbose, accu=acc, return_shift=True
                    )

                    # Guardar shifts de modulación por cuadrante
                    shifts[lambd, q, 0, 0, :] = np.array(srow_mod_q, dtype=float)  # axis=0 (row)
                    shifts[lambd, q, 0, 1, :] = np.array(scol_mod_q, dtype=float)  # axis=1 (col)

                    # # Aplicar shift de modulación SOLO en este cuadrante a ambas cámaras
                    # for cam in range(2):
                    #     for j in range(nmods):
                    #         dy_mod = float(srow_mod_q[j])
                    #         dx_mod = float(scol_mod_q[j])
                    #         aligned[cam, lambd, j] = _apply_patch_shift(
                    #             aligned[cam, lambd, j], y1, y2, x1, x2, dy_mod, dx_mod,
                    #             wrap=True, fill=0
                    #         )

                    srs = ", ".join([f"{float(s):+.3f}" for s in srow_mod_q])
                    scs = ", ".join([f"{float(s):+.3f}" for s in scol_mod_q])
                    print(f"  Q{q}: mod shifts dy=[{srs}]  dx=[{scs}]")

                Q = len(quadrants_roi)
                # Construir shifts por mod para todos los cuadrantes
                shifts_q_per_mod = []
                for j in range(nmods):
                    # para mod j, vector de Q tuples (dy,dx)
                    shifts_for_j = [(float(shifts[lambd, q, 0, 0, j]), float(shifts[lambd, q, 0, 1, j]))
                                    for q in range(Q)]
                    shifts_q_per_mod.append(shifts_for_j)

                # Aplica por cuadrante a CAM1 y CAM2
                for cam in range(2):
                    aligned[cam, lambd] = apply_quadrant_shifts_with_feather_stack(
                        aligned[cam, lambd],             # (4,X,Y)
                        quadrants_roi=quadrants_roi,
                        shifts_q_per_mod=shifts_q_per_mod,
                        margin=10,
                        shift_func=shift_subp,
                        wrap=False, fill=0
                    )


                if debug['activate']:
                    # Visualización comparando pseudo[0] vs pseudo[j] BEFORE/AFTER en un cuadrante
                    # (recalcular pseudo_after aplicando shifts de modulación por cuadrante)
                # Pintar cada mod
                    for j, shifts_mod in enumerate(shifts_q_per_mod, start=1):
                        dy, dx = build_displacement_maps_simple(H, W, quadrants_roi, shifts_mod)
                        plot_maps(dy, dx, title=f"MOD {j}")
                    
                    pseudo_after = pseudo_imgs.copy()
                    for q, (y1, y2, x1, x2) in enumerate(quadrants_roi):
                        dy_vec = shifts[lambd, q, 0, 0, :]  # (4,)
                        dx_vec = shifts[lambd, q, 0, 1, :]
                        # aplica por mod (j) en este cuadrante
                        for j in range(nmods):
                            pseudo_after[j, y1:y2, x1:x2] = shift_subp(
                                pseudo_after[j, y1:y2, x1:x2],
                                shift=[float(dy_vec[j]), float(dx_vec[j])],
                                wrap=False, fill=0
                            )

                    plot_pseudo_mod_alignment(
                        pseudo_before=pseudo_imgs,
                        pseudo_after=pseudo_after,
                        srow_mod=np.zeros(nmods),  # no único; ya que es por cuadrante
                        scol_mod=np.zeros(nmods),
                        zoom_center=debug['zoom_center'],
                        zoom_size=debug['zoom_size'],
                        title_prefix=f"λ={lambd} – pseudo‑modulation alignment (per quadrant)"
                    )

                logging.info("Per‑quadrant camera & modulation alignment applied.")

            if align_sequence == 0:

                # compute shifts for quadrants in camera 0. 
                for q, roises in enumerate(quadrants_roi):
                    y1, y2, x1, x2 = roises
                    patch = data[0, lambd, :, y1:y2, x1:x2]
                    
                    _, srow, scol, _ = realign_subpixel(patch, verbose = verbose, accu = acc, return_shift=True)

                    shifts[lambd, q, 0, 0] = np.array(srow)
                    shifts[lambd, q, 0, 1] = np.array(scol)

                    for npol in range(1,nmods):
                            aligned[0, lambd,npol] = shift_subp(data[0, lambd,npol], shift=[shifts[lambd, q, 0, 0, npol], shifts[lambd, q, 0, 1, npol]], wrap=True, fill=0)

                    # compute shifts for quadrants in camera 0. 
                    for q, roises in enumerate(quadrants_roi):
                        y1, y2, x1, x2 = roises
                        patch = data[1, lambd, :, y1:y2, x1:x2]
                        
                        _, srow, scol, _ = realign_subpixel(patch, verbose = verbose, accu = acc, return_shift=True)

                        shifts[lambd, q, 1, 0] = np.array(srow)
                        shifts[lambd, q, 1, 1] = np.array(scol)

                        for npol in range(1,nmods):
                                aligned[1, lambd,npol] = shift_subp(data[1, lambd,npol], shift=[shifts[lambd, q, 1, 0, npol], shifts[lambd, q, 1, 1, npol]], wrap=True, fill=0)

                    # compute shifts for quadrants from camera 0 (corrected) to camera 1. 
                    for npol in range(nmods):
                        for q, roises in enumerate(quadrants_roi):
                            y1, y2, x1, x2 = roises
                            patch_1 = aligned[0, lambd, npol, y1:y2, x1:x2]
                            patch_2 = aligned[1, lambd, npol, y1:y2, x1:x2] # ojo era rotated

                            _, srow, scol, _ = realign_subpixel(np.array([patch_1,patch_2]), verbose = verbose, accu = acc, return_shift=True)
                            shifts[lambd, q, 1, 0, npol] = srow[1]
                            shifts[lambd, q, 1, 1, npol] = scol[1]

                        # apply shifts for quadrants to camera 1. 

                            aligned[1, lambd,npol] = shift_subp(aligned[1, lambd,npol], shift=[shifts[lambd, q, 1, 0, npol], shifts[lambd, q, 1, 1, npol]], wrap=True, fill=0) # ojo era rotated

            else:
                pass

        else:
            raise ValueError("Quadrants parameter not recognized. It must be >= 0.")

    tac = time.time()

    logging.info(f"Alignment finished in {round(tac - tic, 3)} s.")

    if returnshifts:
        if onelambda:
            return aligned[:, 0], shifts[0]
        else:    
            return aligned, shifts, err
    else:
        if onelambda:
            return aligned[:, 0]
        else:    
            return aligned

from scipy.signal import fftconvolve
from scipy.ndimage import center_of_mass, gaussian_filter
import numpy as np

def moffat_psf(size, alpha, beta):
    ax = np.arange(-(size // 2), size // 2 + 1)
    xx, yy = np.meshgrid(ax, ax)
    r2 = xx**2 + yy**2

    psf = (1 + (r2 / alpha**2))**(-beta)
    psf /= psf.sum()
    return psf

def azimuthal_average(image, center):
    y, x = np.indices(image.shape)
    r = np.sqrt((x - center[1])**2 + (y - center[0])**2)

    r_int = r.astype(int)

    tbin = np.bincount(r_int.ravel(), image.ravel())
    nr = np.bincount(r_int.ravel())

    radial_profile = tbin / np.maximum(nr, 1)

    radius = np.arange(len(radial_profile))

    return radial_profile, radius

def detect_disk(image, smooth_sigma=2):
    p_low = np.percentile(image, 5)
    p_high = np.percentile(image, 95)

    background = p_low
    image_corr = image - background
    image_corr[image_corr < 0] = 0
    if p_high > p_low:
        image_corr /= (p_high - p_low)
    image_smooth = gaussian_filter(image_corr, sigma=smooth_sigma)
    cy, cx = center_of_mass(image_smooth)
    y, x = np.indices(image.shape)
    r = np.sqrt((x - cx)**2 + (y - cy)**2)
    gy, gx = np.gradient(image_smooth)
    grad = np.sqrt(gx**2 + gy**2)
    prof_grad, r_arr = azimuthal_average(grad, (cy, cx))

    R = r_arr[np.argmax(prof_grad)]

    i_peak = np.argmax(prof_grad)

    if 1 < i_peak < len(prof_grad) - 2:
        y0 = prof_grad[i_peak-1:i_peak+2]
        x0 = r_arr[i_peak-1:i_peak+2]

        # ajuste parabólico simple
        try:
            coeff = np.polyfit(x0, y0, 2)
            R = -coeff[1] / (2 * coeff[0])
        except:
            pass

    return (cy, cx, R)

def estimate_epsilon(image, target_contrast=0.15, mask_radius_factor=0.3, epsilon = 0):
    cy, cx, radius = detect_disk(image)

    y, x = np.indices(image.shape)
    r = np.sqrt((x - cx)**2 + (y - cy)**2)

    mask = r <= (mask_radius_factor * radius)

    Iq = image[mask]
    contrast = np.std(Iq) / np.mean(Iq)

    if epsilon == 0:
        eps = max(0, 1 - contrast / target_contrast)
    else:
        eps = epsilon

    return np.clip(eps, 0.01, 0.4), mask

def correct_straylight(data,
                    psf_size=31,
                    alpha=4.0,
                    beta=2.5,
                    n_iter=2,
                    target_contrast=0.15,
                    epsilon = 0.10,
                    mask_radius_factor=0.3,
                    wave_ref=0):

    """
    Wavelength-independent stray light correction

    Parameters
    ----------
    alpha, beta : PSF parameters (fixed)
    epsilon     : estimated once from a reference wavelength if epsilon == 0. If not, just take the value.

    alpha ≈ 3 – 6 pixels
    beta  ≈ 2 – 3
    epsilon ≈ 0.05 – 0.20

    """

    nw, ns, nx, ny = data.shape

    # --- Estimate epsilon ONLY ONCE ---
    epsilon, mask = estimate_epsilon(
        data[wave_ref, 0],
        target_contrast=target_contrast,
        epsilon = epsilon,
        mask_radius_factor=mask_radius_factor
    )

    # --- Single PSF ---
    psf = moffat_psf(psf_size, alpha, beta)

    corrected = data.copy()

    for _ in range(n_iter):
        for w in range(nw):
            for s in range(ns):
                img = corrected[w, s]

                blurred = fftconvolve(img, psf, mode='same')

                corrected[w, s] = img - epsilon * (blurred - img)

    return corrected, psf, mask


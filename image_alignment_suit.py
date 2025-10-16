import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import affine_transform
from scipy.optimize import minimize

def extract_patch(image, center, size=20):
    y, x = np.round(center).astype(int)
    half = size // 2
    return image[y-half:y+half, x-half:x+half]

def cost_function_pixelwise(params, image0, image1, positions, size=40):
    angle_deg, t_y, t_x, c_y, c_x, scale_x, scale_y, shear_x, shear_y = params
    t = np.array([t_y, t_x])
    center = np.array([c_y, c_x])
    transformed_image1 = apply_transform(image1, angle_deg, t, center,
                                         scale_x, scale_y, shear_x, shear_y)

    rms_list = []
    for p0 in positions:
        patch0 = extract_patch(image0, p0, size=size)
        patch1 = extract_patch(transformed_image1, p0, size=size)
        if patch0.shape == patch1.shape:
            rms = np.std((patch0 - patch1)**2)
            rms_list.append(rms)
    return np.mean(rms_list)

def apply_transform(image, angle_deg, t, center,
                    scale_x=1.0, scale_y=1.0, shear_x=0.0, shear_y=0.0):
    angle_rad = np.radians(angle_deg)
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
    R = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
    Distortion = np.array([[scale_x, shear_x], [shear_y, scale_y]])
    A = Distortion @ R
    A_inv = np.linalg.inv(A)
    offset = center - A_inv @ (center + t)

    return affine_transform(
        image, matrix=A_inv, offset=offset, order=3, mode='nearest', cval=0.0
    )

def normalize_image(img):
    p1, p99 = np.percentile(img, (1, 99))
    img_clipped = np.clip(img, p1, p99)
    return (img_clipped - img_clipped.min()) / (img_clipped.max() - img_clipped.min())


def image_alignment_affine(image0,image1, init_params = [0.5, 2., 2.,
                   1000., 1000., 1.0, 1.0, 0.0, 0.0]):

    init_center = image1.shape

    bounds = [
        (- 0.0, + 0.0),
        ( - 1,  + 1),
        ( - 1,  + 1),
        (init_center[0]//2 - 200, init_center[0]//2 + 200),
        (init_center[1]//2 - 200, init_center[1]//2 + 200),
        (0.995, 1.005), (0.995, 1.005),
        (-0.005, 0.005), (-0.005, 0.005)
    ]
    # Create a 3x3 grid of positions distributed along the image
    h, w = image0.shape
    ys = np.linspace(h * 0.2, h * 0.8, 3) - 100
    xs = np.linspace(w * 0.2, w * 0.8, 3) + 100
    positions = np.array([[y, x] for y in ys for x in xs])

    res = minimize(cost_function_pixelwise, init_params,
                   args=(image0, image1, positions), bounds=bounds)

    (opt_angle, opt_ty, opt_tx, opt_cy, opt_cx,
     opt_sx, opt_sy, opt_shx, opt_shy) = res.x

    print("\nOptimized transformation parameters:")
    print(f"  Rotation angle (deg): {opt_angle:.4f}")
    print(f"  Translation (t_y, t_x): ({opt_ty:.2f}, {opt_tx:.2f})")
    print(f"  Center of rotation (c_y, c_x): ({opt_cy:.2f}, {opt_cx:.2f})")
    print(f"  Scale (x, y): ({opt_sx:.4f}, {opt_sy:.4f})")
    print(f"  Shear (x, y): ({opt_shx:.5f}, {opt_shy:.5f})")
    print(f"  Minimal RMS error: {res.fun:.4f}")

    aligned_image1 = apply_transform(
        image1, opt_angle, np.array([opt_ty, opt_tx]),
        np.array([opt_cy, opt_cx]),
        scale_x=opt_sx, scale_y=opt_sy, shear_x=opt_shx, shear_y=opt_shy
    )

    plt.figure(figsize=(10, 10))
    plt.subplot(1, 1, 1)
    plt.title("Difference (normalized)")
    diff_img = np.abs(image0 - aligned_image1)
    diff_norm = normalize_image(diff_img)
    plt.imshow(diff_img, cmap='inferno')
    plt.colorbar()
    plt.show()

    return aligned_image1


import numpy as np
from alignment import *

def reshape_into_16_quadrants(images, nlambda, nmods):
    """
    Function to reshape an observing mode into 16 quadrant along the FoV (4 x 4):
    Inputs: 
        - images (np.array) : Array contaning the obs mode. (Ncams x Nlambda x Nmods x 1416 x 1416)
        - nlambda (int): Number of wavelengths
        - nmods (int) : Number of modulations. 
    Outputs:
        - reshaped (np.array) : Array with an additional axis to go over the quadrants (Ncams x Nlambda x Nmods x 16 x 354 x 35)
    """
    # Reshape the last two dimensions into a 4x4 grid of (354, 354)
    reshaped = images.reshape(2, nlambda, nmods, 4, 354, 4, 354)
    # Rearrange axes to group quadrants into a single dimension
    return reshaped.transpose(0, 1, 2, 3, 5, 4, 6).reshape(2, nlambda, nmods, 16, 354, 354)

def align_quadrants(data, acc = 0.01, verbose = False):
    """
    Function to align the quadrants separately
    Inputs: 
        - data (np.array) : Array contaning the obs mode. (Ncams x Nlambda x Nmods x 1416 x 1416)
        - acc (int, default : 0.01): accuracy for the alignment routine
        - verbose (Boolean, default : False) : Print info on terminal 
    Outputs:
        - aligned (np.array) : Array with an additional axis to go over the quadrants (Ncams x Nlambda x Nmods x 16 x 354 x 35)
        - shifts (list) : shifts performed to each camera, modulation and wavelength
    """
    shape = np.shape(data)
    nlambda = shape[1]
    nmods = shape[2]
    nquads = shape[3]

    shifts = np.zeros((nlambda, 2, 2, nmods, nquads))

    aligned = np.zeros(np.shape(data))

    for lambd in range(nlambda):

        print(f"\nAligning wavelengh: {lambd}/{nlambda}")
        print(f"-------------------------------------")

        for quad in range(nquads):

            print(f"\nProcessing Q{quad}... \nModulations of cam 1 alignment...")
            mods_aligned, srow, scol = realign_subpixel(data[0, lambd, :, quad], verbose = verbose, accu = acc, return_shift=True)

            shifts[lambd, 0, 0, :, quad] = srow
            shifts[lambd, 0, 1, :, quad] = scol

            aligned[0, lambd, :, quad] = mods_aligned

            for mod in range(nmods):
                print(f"Cam 2 of M{mod}...")
                cams_aligned, srow, scol = realign_subpixel(np.array([mods_aligned[mod], data[1, lambd, mod, quad]]), verbose = verbose, accu = acc, return_shift=True )

                shifts[lambd, 1, 0, mod, quad] = srow[1]
                shifts[lambd, 1, 1, mod, quad] = scol[1]

                aligned[1, lambd, mod, quad] = cams_aligned[1]

    return aligned, shifts

def _apply_patch_shift(img, y1, y2, x1, x2, dy, dx, wrap=False, fill=0):
    """
    Desplaza SOLO el parche [y1:y2, x1:x2] de una imagen 2D y lo pega de vuelta.
    """
    patch = img[y1:y2, x1:x2]
    patch_shifted = shift_subp(patch, shift=[dy, dx], wrap=wrap, fill=fill)
    out = img.copy()
    out[y1:y2, x1:x2] = patch_shifted
    return out

def _apply_patch_shift_stack(stack4, y1, y2, x1, x2, dy, dx, wrap=False, fill=0):
    """
    Igual que _apply_patch_shift pero para un stack (4, H, W) de modulaciones.
    Aplica el mismo (dy, dx) a cada capa.
    """
    out = stack4.copy()
    for j in range(stack4.shape[0]):
        patch = stack4[j, y1:y2, x1:x2]
        out[j, y1:y2, x1:x2] = shift_subp(patch, shift=[dy, dx], wrap=wrap, fill=fill)
    return out

def rotate_camera2(data, theta = 0.065, onelambda = False):
    """
    Function to rotate the camera 2 from an obs mode

    Inputs:
        - data (np.array) : Array contaning the obs mode. (Ncams x Nlambda x Nmods x Nx x Ny)
        - theta (float, default : 0.065): Angle of rotation
        - onelambda (Boolean, default : False): Set to true if only one lambda is used (array of shape Ncam x Nmod x Nx x Ny) 
    Outputs:
        - Rotation (np.array) : Same array as data with cam2 rotated by theta 
    """

    if onelambda:
        data = data[:, np.newaxis] # To allow for only one lamdba.

    # Get shape for data
    shape = np.shape(data)
    nlambda = shape[1]
    nmods = shape[2]

    print("Computing camera 2 rotation...")

    rotated = np.copy(data)
    for lambd in range(nlambda):
        for mod in range(nmods):    
            rotated[1, lambd, mod] = rotate(data[1, lambd, mod], theta, reshape=False, order = 2)  
    print("Rotation finished.\n")  

    return rotated

def filter_and_rotate(data, theta = 0.0655, verbose = False, filterflag = True, zkes = np.zeros(21)):

    """
    Function to filter an obs mode and rotate camera 2

    Inputs:
        - data (np.array) : Array contaning the obs mode. (Ncams x Nlambda x Nmods x Nx x Ny)
        - theta (float, default : 0.0655): Angle of rotation
        - verbose (Boolean, ddefault : False) : Print info on terminal.
        - filterflag (Boolean, default : True) : Set to False to skip Fourier filtration 
        - zkes (np.array, default : np.zeros(21)): Zernike's array to use for the filtration. 
    Outputs:
        - Filtered and Rotated (np.array) : Same array as data filtrated and with cam2 rotated  
    """

    # Get shape for data
    shape = np.shape(data)
    nlambda = shape[1]
    nmods = shape[2]

    filtered_n_rotated  = np.zeros(shape)

    if filterflag and verbose:
        print("Noise filtration and camera 2 rotation...")
    elif verbose:
        print("Cam 2 rotation...")

    for lambd in range(nlambda):
        print(f"Procesing wavelength: {lambd + 1} / {nlambda}")
        
        for mod in range(nmods):
            # Apply noise filter

            if filterflag:
                filtered_n_rotated[0, lambd, mod], _ = restore_ima(data[0, lambd, mod], zkes)
                cam2_filtered, _ = restore_ima(data[1, lambd, mod], zkes)
                filtered_n_rotated[1, lambd, mod] = rotate(cam2_filtered, theta, reshape=False, order = 2)
            else:
                filtered_n_rotated[0, lambd, mod] = data[0, lambd, mod]
                filtered_n_rotated[1, lambd, mod] = rotate(data[1, lambd, mod], theta, reshape=False, order = 2)
    
    return filtered_n_rotated


def demodulate_quadrants(data, nlambda, nmods, filt, nquads = 16, Np_quad = 354):
    """
    Function to perform the demodulation of the observation mode separated in quadrants. 
    Inputs: 
        - data (np.array) : Array contaning the obs mode. (Ncams x Nlambda x Nmods x Nx x Ny)
        - nlambda (int): Number of wavelengths
        - nmods (int) : Number of modulations
        - filt (str) : Filter to demodulate (517, 525.02 or 525.06)
        - nquads (int, default : 16) : Number of quadrants
        - Np_quad (int, defaulr : 354) : Pixel soize of quadrant
    Outputs:
        - dual_beamed (np.array) : Demodulated data with cameras combined (Nlambda x Nmods x Nx x Ny).
        - demodulated (np.array) : Demodulated data with cameras not yet combined (Ncams x Nlambda x Nmods x Nx x Ny). 
    """

    demod = np.zeros((2, nlambda, nmods, nquads, Np_quad, Np_quad))
    dual = np.zeros((nlambda, nmods, nquads, Np_quad, Np_quad))

    for quad in range(nquads):

        du, dem = demodulate(data[:, :, :, quad], nmods, nlambda, filt)
        demod[:, :, :, quad] = dem
        dual[:, :, quad] = du

    return dual, demod
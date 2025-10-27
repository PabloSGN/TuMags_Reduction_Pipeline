# ---------------------------- DESCRIPTION --------------------------------------- #
"""

Module with all the functions related to the alignment of the observation modes. 

Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ IMPORTS ----------------------------------------- #

import numpy as np
import time
from matplotlib import pyplot as plt
from scipy.fftpack import fftshift, ifftshift, fft2, ifft2
from scipy.ndimage import rotate

# Own functions
from pd_functions_v22 import restore_ima
from image_filtering import filter_frecuencies

# ------------------------------  AUX FUNCTIONS  --------------------------------- # 

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

def align_obsmode(data, acc = 0.01, verbose = False, theta = 0.0655, filterflag = False, 
                  onelambda = False, returnshifts = True,roi = [0,-1,0,-1], quadrants = 0):
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
    # rotated = np.copy(data)

    if quadrants != 0:
        H, W = data.shape[-1], data.shape[-1]   # dimensiones de la imagen

        quadrants_roi = _generar_cuadrantes_nx_n(H, W, quadrants)

        # for idx, (y1, y2, x1, x2) in enumerate(quadrants_roi):
        #     print(f"ROI {idx}: y={y1}:{y2}, x={x1}:{x2}")
        
        shifts = np.zeros((nlambda, quadrants*quadrants, 2, 2, nmods), dtype=float)
    else:
        shifts = np.zeros((nlambda, 2, 2, nmods))

    err = []
    for lambd in range(nlambda):
        print(f"Aligning wavelength: {lambd + 1}/{nlambda}")
        if quadrants == 0:
            if verbose:
                print(f"Aligning wavelengh: {lambd + 1}/{nlambda}")
                print(f"Shifts for cam 1 - modulation alignment")

            _, srow, scol, error = realign_subpixel(data[0, lambd,:,roi[0]:roi[1],roi[2]:roi[3]], verbose = verbose, accu = acc, return_shift=True)
            err.append(error)

            shifts[lambd, 0, 0] = srow
            shifts[lambd, 0, 1] = scol

            for nm in range(nmods-1):
                aligned[0, lambd,nm + 1] = shift_subp(data[0, lambd,nm+1], shift=[srow[nm + 1], scol[nm + 1]], wrap=True, fill=0)

            if verbose:
                print("Shifts of camera 2 alignment")
            for mod in range(nmods):
                if verbose:
                    print(f"mod -> {mod}...")

                _, srow, scol,error = realign_subpixel(np.array([aligned[0,lambd,mod,roi[0]:roi[1],roi[2]:roi[3]], data[1, lambd, mod,roi[0]:roi[1],roi[2]:roi[3]]]), verbose = verbose, accu = acc, return_shift=True)
                err.append(error)

                shifts[lambd, 1, 0, mod] = srow[1]
                shifts[lambd, 1, 1, mod] = scol[1]

                aligned[1, lambd,mod] = shift_subp(data[1, lambd,mod], shift=[srow[1], scol[1]], wrap=True, fill=0)
        else:

            both_cams = 2
            # compute shifts for quadrants in camera 0. 
            for q, roises in enumerate(quadrants_roi):
                y1, y2, x1, x2 = roises
                patch = data[0, lambd, :, y1:y2, x1:x2]
                
                _, srow, scol, _ = realign_subpixel(patch, verbose = verbose, accu = acc, return_shift=True)

                shifts[lambd, q, 0, 0] = np.array(srow)
                shifts[lambd, q, 0, 1] = np.array(scol)

                for npol in range(1,nmods):
                        aligned[0, lambd,npol] = shift_subp(data[0, lambd,npol], shift=[shifts[lambd, q, 0, 0, npol], shifts[lambd, q, 0, 1, npol]], wrap=True, fill=0)

            if both_cams == 1 or both_cams == 2:
                # compute shifts for quadrants in camera 0. 
                for q, roises in enumerate(quadrants_roi):
                    y1, y2, x1, x2 = roises
                    patch = data[1, lambd, :, y1:y2, x1:x2]
                    
                    _, srow, scol, _ = realign_subpixel(patch, verbose = verbose, accu = acc, return_shift=True)

                    shifts[lambd, q, 1, 0] = np.array(srow)
                    shifts[lambd, q, 1, 1] = np.array(scol)

                    for npol in range(1,nmods):
                            aligned[1, lambd,npol] = shift_subp(data[1, lambd,npol], shift=[shifts[lambd, q, 1, 0, npol], shifts[lambd, q, 1, 1, npol]], wrap=True, fill=0)

            if both_cams == 0 or both_cams == 2:

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
            elif both_cams == 1:
                # compute shifts for quadrants from camera 0 (corrected) to camera 1. 
                for npol in range(nmods):
                    patch_1 = aligned[0, lambd, npol, roi[0]:roi[1],roi[2]:roi[3]]
                    patch_2 = aligned[1, lambd, npol, roi[0]:roi[1],roi[2]:roi[3]]

                    _, srow, scol, _ = realign_subpixel(np.array([patch_1,patch_2]), verbose = verbose, accu = acc, return_shift=True)

                    aligned[1, lambd,npol] = shift_subp(aligned[1, lambd,npol], shift=[srow[1], scol[1]], wrap=True, fill=0)
            else:
                pass

    # err = []
    # for lambd in range(nlambda):
    #     print(f"Aligning wavelength: {lambd + 1}/{nlambda}")
    #     if quadrants == 0:
    #         if verbose:
    #             print(f"Aligning wavelengh: {lambd + 1}/{nlambda}")
    #             print(f"Shifts for cam 1 - modulation alignment")

    #         _, srow, scol, error = realign_subpixel(rotated[0, lambd,:,roi[0]:roi[1],roi[2]:roi[3]], verbose = verbose, accu = acc, return_shift=True)
    #         err.append(error)

    #         shifts[lambd, 0, 0] = srow
    #         shifts[lambd, 0, 1] = scol

    #         for nm in range(nmods-1):
    #             aligned[0, lambd,nm + 1] = shift_subp(rotated[0, lambd,nm+1], shift=[srow[nm + 1], scol[nm + 1]], wrap=True, fill=0)

    #         if verbose:
    #             print("Shifts of camera 2 alignment")
    #         for mod in range(nmods):
    #             if verbose:
    #                 print(f"mod -> {mod}...")

    #             _, srow, scol,error = realign_subpixel(np.array([aligned[0,lambd,mod,roi[0]:roi[1],roi[2]:roi[3]], rotated[1, lambd, mod,roi[0]:roi[1],roi[2]:roi[3]]]), verbose = verbose, accu = acc, return_shift=True)
    #             err.append(error)

    #             shifts[lambd, 1, 0, mod] = srow[1]
    #             shifts[lambd, 1, 1, mod] = scol[1]

    #             aligned[1, lambd,mod] = shift_subp(rotated[1, lambd,mod], shift=[srow[1], scol[1]], wrap=True, fill=0)
    #     else:

    #         both_cams = 2
    #         # compute shifts for quadrants in camera 0. 
    #         for q, roises in enumerate(quadrants_roi):
    #             y1, y2, x1, x2 = roises
    #             patch = rotated[0, lambd, :, y1:y2, x1:x2]
                
    #             _, srow, scol, _ = realign_subpixel(patch, verbose = verbose, accu = acc, return_shift=True)

    #             shifts[lambd, q, 0, 0] = np.array(srow)
    #             shifts[lambd, q, 0, 1] = np.array(scol)

    #             for npol in range(1,nmods):
    #                     aligned[0, lambd,npol] = shift_subp(rotated[0, lambd,npol], shift=[shifts[lambd, q, 0, 0, npol], shifts[lambd, q, 0, 1, npol]], wrap=True, fill=0)

    #         if both_cams == 1 or both_cams == 2:
    #             # compute shifts for quadrants in camera 0. 
    #             for q, roises in enumerate(quadrants_roi):
    #                 y1, y2, x1, x2 = roises
    #                 patch = rotated[1, lambd, :, y1:y2, x1:x2]
                    
    #                 _, srow, scol, _ = realign_subpixel(patch, verbose = verbose, accu = acc, return_shift=True)

    #                 shifts[lambd, q, 1, 0] = np.array(srow)
    #                 shifts[lambd, q, 1, 1] = np.array(scol)

    #                 for npol in range(1,nmods):
    #                         aligned[1, lambd,npol] = shift_subp(rotated[1, lambd,npol], shift=[shifts[lambd, q, 1, 0, npol], shifts[lambd, q, 1, 1, npol]], wrap=True, fill=0)

    #         if both_cams == 0 or both_cams == 2:

    #             # compute shifts for quadrants from camera 0 (corrected) to camera 1. 
    #             for npol in range(nmods):
    #                 for q, roises in enumerate(quadrants_roi):
    #                     y1, y2, x1, x2 = roises
    #                     patch_1 = aligned[0, lambd, npol, y1:y2, x1:x2]
    #                     patch_2 = aligned[1, lambd, npol, y1:y2, x1:x2] # ojo era rotated

    #                     _, srow, scol, _ = realign_subpixel(np.array([patch_1,patch_2]), verbose = verbose, accu = acc, return_shift=True)
    #                     shifts[lambd, q, 1, 0, npol] = srow[1]
    #                     shifts[lambd, q, 1, 1, npol] = scol[1]

    #                 # apply shifts for quadrants to camera 1. 

    #                     aligned[1, lambd,npol] = shift_subp(aligned[1, lambd,npol], shift=[shifts[lambd, q, 1, 0, npol], shifts[lambd, q, 1, 1, npol]], wrap=True, fill=0) # ojo era rotated
    #         elif both_cams == 1:
    #             # compute shifts for quadrants from camera 0 (corrected) to camera 1. 
    #             for npol in range(nmods):
    #                 patch_1 = aligned[0, lambd, npol, roi[0]:roi[1],roi[2]:roi[3]]
    #                 patch_2 = aligned[1, lambd, npol, roi[0]:roi[1],roi[2]:roi[3]]

    #                 _, srow, scol, _ = realign_subpixel(np.array([patch_1,patch_2]), verbose = verbose, accu = acc, return_shift=True)

    #                 aligned[1, lambd,npol] = shift_subp(aligned[1, lambd,npol], shift=[srow[1], scol[1]], wrap=True, fill=0)
    #         else:
    #             pass

    tac = time.time()

    if verbose:
        print(f"Alignment finished in {round(tac - tic, 3)} s.")

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
        

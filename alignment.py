# ---------------------------- DESCRIPTION --------------------------------------- #

"""
author: Pablo Santamarina Guerrero(pablosantamarinag@gmail.com) 
Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ IMPORTS ----------------------------------------- #

import numpy as np
import time
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from scipy.fftpack import fftshift, ifftshift, fft2, ifft2
from scipy.ndimage import rotate

# Own functions
from pd_functions_v22 import restore_ima
from image_filtering import filter_frecuencies

# Config
margin = 10 # Margin from fieldstop

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
    Nr=np.fft.ifftshift(np.arange(-np.fix(nr/2),np.ceil(nr/2)))
    Nc=np.fft.ifftshift(np.arange(-np.fix(nc/2),np.ceil(nc/2)))
    CC=np.fft.ifft2(FTpad(F*np.conj(G),2*nr))
    CCabs=np.abs(CC)
    ind = np.unravel_index(np.argmax(CCabs, axis=None), CCabs.shape)
    CCmax=CC[ind]*nr*nc
    Nr2=np.fft.ifftshift(np.arange(-np.fix(nr),np.ceil(nr)))
    Nc2=np.fft.ifftshift(np.arange(-np.fix(nc),np.ceil(nc)))
    row_shift=Nr2[ind[0]]/2
    col_shift=Nr2[ind[1]]/2

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
    np.fft.ifftshift(np.arange(0,nc).T-np.floor(nc/2)),np.arange(0,n_out)-coff))

    kernr=np.exp((-1j*2*np.pi/(nr*kappa))*np.outer(\
    np.arange(0,n_out)-roff,np.fft.ifftshift(np.arange(0,nr).T-np.floor(nr/2))))
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
    IM=np.fft.fftshift(IM)
    IMout=np.pad(IM,((pd,pd),(pd,pd)),'constant')
    IMout=np.fft.ifftshift(IMout)*Nout*Nout/(Nin*Nin)
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
        return ima_aligned, row_shifts, col_shifts
    else:
        return ima_aligned

def find_fieldstop(cam1 = None, verbose = False, plot_flag = False, margin = margin):

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

def align_obsmode(data, acc = 0.001, verbose = False, theta = 0.0655, filterflag = True, onelambda = False, zkes = np.zeros(21)):

    if onelambda:
        data = data[:, np.newaxis] # To allow for only one lamdba.

    shape = np.shape(data)
    nlambda = shape[1]
    nmods = shape[2]

    shifts = np.zeros((nlambda, 2, 2, nmods))
    aligned = np.zeros(np.shape(data))

    if filterflag:
        filtered = filter_frecuencies(data)
        rotated = rotate_camera2(filtered, theta = theta)
    else:
        rotated = rotate_camera2(data, theta = theta)

    # Start by applying noise and rotation
    """filtered_and_rotated = filter_and_rotate(data, verbose = verbose, 
                                             theta=theta, filterflag = filterflag, zkes= zkes)""" 

    for lambd in range(nlambda):

        print(f"Aligning wavelengh: {lambd}/{nlambda}")

        print(f"Modulations of cam 1 alignment...")
        mods_aligned, srow, scol = realign_subpixel(rotated[0, lambd], verbose = verbose, accu = acc, return_shift=True)

        shifts[lambd, 0, 0] = srow
        shifts[lambd, 0, 1] = scol

        aligned[0, lambd] = mods_aligned

        print("Aligning cam2...")
        for mod in range(nmods):
            print(f"mod -> {mod}...")
            cams_aligned, srow, scol = realign_subpixel(np.array([mods_aligned[mod], rotated[1, lambd, mod]]), verbose = verbose, accu = acc, return_shift=True )

            shifts[lambd, 1, 0, mod] = srow[1]
            shifts[lambd, 1, 1, mod] = scol[1]

            aligned[1, lambd, mod] = cams_aligned[1]

    if onelambda:
        return aligned[:, 0], shifts[0]
    else:    
        return aligned, shifts

def reshape_into_16_quadrants(images, nlambda, nmods):
    # Reshape the last two dimensions into a 4x4 grid of (354, 354)
    reshaped = images.reshape(2, nlambda, nmods, 4, 354, 4, 354)
    # Rearrange axes to group quadrants into a single dimension
    return reshaped.transpose(0, 1, 2, 3, 5, 4, 6).reshape(2, nlambda, nmods, 16, 354, 354)

def align_quadrants(data, acc = 0.01, verbose = False):

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
        
"""
This module contains the functions to correct for crosstalk
through the method described in Jaeggli et al. 2022. 
(https://doi.org/10.3847/1538-4357/ac6506)

https://github.com/sajaeggli/adhoc_xtalk/blob/main/Jaeggli_etal_2022ApJ_AdHoc_Xtalk.ipynb

"""
import numpy as np
from matplotlib import colors, pyplot as plt
from glob import glob
from scipy.optimize import minimize
import pandas as pd

# These are the model and minimization functions as defined in the paper

# Functions for the diattenuation modeling
def polmodel1(D,theta,chi):
    dH = D*np.cos(chi)*np.sin(theta)
    d45 = D*np.sin(chi)*np.sin(theta)
    dR = D*np.cos(theta)
    A = np.sqrt(1. - dH**2 - d45**2 - dR**2)
    
    mat1 = np.array([
        [ 1., dH, d45, dR], 
        [ dH,  A,  0., 0.], 
        [d45, 0.,   A, 0.],
        [ dR, 0.,  0.,  A]], dtype='double')
    
    mat2 = np.array([
        [0.,     0.,     0.,     0.],
        [0.,  dH**2, d45*dH,  dH*dR],
        [0., d45*dH, d45**2, d45*dR],
        [0.,  dH*dR, d45*dR,  dR**2]], dtype='double')
    
    return( mat1 + (1-A)/D**2*mat2 )

#Function that returns the merit function from the Stokes vector
#and the values of the parameters that define the diattenuation matrix
def fitfunc1(param, stokesin):
    D = param[0]
    theta = param[1]
    chi = param[2]
    
    # Keep diattenuation value in range
    if D>=1:
        D=0.999999
        
    if D<=-1:
        D = -0.999999
    
    #Computes Mueller Matrix (MM) and its inverse (iMM)
    MM = polmodel1(D, theta, chi)
    iMM = np.linalg.inv(MM)


    #Computes the merit function
    out = minimize_for_model1(iMM,stokesin)

    return(out)

#Function that computes the merit function for a given Mueller matrix
def minimize_for_model1(iMM,bs):
    # apply a mueller matrix (rotation) to a 3D stokes vector 
    # (npixel,wavelength,4)
    new_stokes = np.einsum('ij,abj->abi',iMM, np.squeeze(bs))

    
    # Minimization criteria
    out = np.abs(np.sum(new_stokes[:,:,0]*new_stokes[:,:,3],axis=1)) + \
          np.abs(np.sum(new_stokes[:,:,0]*new_stokes[:,:,2],axis=1)) + \
          np.abs(np.sum(new_stokes[:,:,0]*new_stokes[:,:,1],axis=1))

    # sum over spatial positions
    out = np.sum(out)
    
    return(out)


# Function for the retarder modeling
def polmodel2(theta, delta):
    St = np.sin(theta)
    Ct = np.cos(theta)
    Sd = np.sin(delta)
    Cd = np.cos(delta)
 
    MM1 = np.array([
        [1.,  0., 0., 0.],
        [0.,  Ct, St, 0.],
        [0., -St, Ct, 0.],
        [0.,  0., 0., 1.]
    ], dtype='double')
    
    MM2 = np.array([
        [1., 0.,  0., 0.],
        [0., 1.,  0., 0.],
        [0., 0.,  Cd, Sd],
        [0., 0., -Sd, Cd]
    ], dtype='double')
    
    MM = np.einsum('ij,jk', MM1, MM2)
    return(MM)

#Function that builds the merit function from the Stokes vector
# and the values of the parameters that define the retarder matrix    
def fitfunc2(fitangles, stokesin):
    theta = fitangles[0]
    delta = fitangles[1]
    
    MM = polmodel2(theta, delta)
    iMM = np.linalg.inv(MM)

    out = minimize_for_model2(iMM, stokesin)

    return(out)

#Function that computes the merit function for a given Mueller matrix
def minimize_for_model2(iMM,bs):
    new_stokes = np.einsum('ij,abj->abi',iMM, np.squeeze(bs))
    
    # Minimization criteria
    out = np.sum(new_stokes[:,:,3],axis=1)**2 +\
          np.abs(np.sum(new_stokes[:,:,3]*new_stokes[:,:,2],axis=1)) +\
          np.abs(np.sum(new_stokes[:,:,3]*new_stokes[:,:,1],axis=1))
    
    # sum over spatial positions
    out = np.sum(out)
    
    return(out)

def fit_mueller_matrix(data,pthresh=0.02,norm=False,
                       region=[200,1200,200,1200],
                       last_wvl=None,plots=False):
    """
    This function fits the best diattenuation Mueller matrix that
    minimizes the correlation between Stokes I to Q, U and V along
    the spectral line of interest over a weakly polarized region.
    Input:
        data: 4D array with the Stokes parameters after dual beam
            with dimensions (wavelength, Stokes, x, y)
        pthresh: Stokes V polarization threshold to determine weak/strong
          polarization regions. Default: 0.02
        norm: True or False. Normalization of Stokes components to
          median of Stokes I at first wavelength.
        region: List of the type (x0,xf,y0,yf) with region of the
          image to be considered for the fit. Normally, a central
          region not affected by edge artifacts arising after 
          alignment/rotation. Default: [200,1200,200,1200].
        last_wvl: -1 (Mg I) or None. Last wavelength to be considered
          to compute the merit function. Default: None
        plots: True or false. If True, plots the fractional 
            polarization map 
    Output:
        datarest: 4D array with the corrected Stokes parameters
        MM1a: Diattenuation Mueller matrix that converts the
            "real" Stokes parameters into the "observed" ones.
    """
    #Reorder axis to convert into dimensions: [x,y,wavelength,stokes]
    data=np.moveaxis(data,0,-1)
    data=np.moveaxis(data,0,-1)


    if norm is True:
        #Normalization of data
        norm_factor=np.median(data[:,:,0,0])
        data=data/norm_factor

    ## Crop data to avoid edge effects arising from alignment/rotation
    full_data=data.copy()
    data=data[region[0]:region[1],region[2]:region[3],:,:]

    #Last wavelength to be considered in the minimization
    data=data[:,:,:last_wvl,:]

    # Choose initial guess parameters for the diattenuation minimization
    D = 0.5
    theta = 0.
    chi = 0.
    initial_guess = (D, theta, chi)

    #Apply threshold to V map corrected from offset.
    V_mean_corr=data[:,:,:,3]-np.mean(data[:,:,0,3],axis=(0,1))
    pmap = np.max(np.abs(V_mean_corr)/data[:,:,:,0], axis=2)
    notpolar = np.argwhere(pmap < pthresh)
    nyidx = notpolar[:,0]
    nzidx = notpolar[:,1]

    # Use just the region with weak polarization
    weak_region = data[nyidx,nzidx,:,:] #do selection for only strong polarization signals

    if plots is True:
        #Plot fractional polarization map
        fig,ax=plt.subplots(figsize=(8,8))
        plot=ax.imshow(pmap)
        ax.set_title('Fractional polarization map')
        plt.colorbar(plot)

        #Plot original data at wavelength 0 and contour of weak/stron regions
        fig,axs=plt.subplots(2,2,layout='constrained',figsize=(10,10))
        axs[0,0].imshow(data[:,:,0,0],cmap='gray')
        axs[0,0].set_title('Stokes I')
        axs[0,1].imshow(data[:,:,0,1],cmap='gray')
        axs[0,1].set_title('Stokes Q')
        axs[1,0].imshow(data[:,:,0,2],cmap='gray')
        axs[1,0].set_title('Stokes U')
        axs[1,1].imshow(data[:,:,0,3],cmap='gray')
        axs[1,1].set_title('Stokes V')
        for i in range(2):
            for j in range(2):
                axs[i,j].contour(pmap, [pthresh], colors='green', 
                                linewidths=0.75)
        plt.show()
        plt.close()        
    #Minimize merit function
    result = minimize(fitfunc1, initial_guess, args=weak_region)


    # Apply correction for I<->QUV cross-talk
    MM1a = polmodel1(result.x[0],result.x[1], result.x[2])
    iMM1a = np.linalg.inv(MM1a)
    data_corrected =  np.einsum('ij,abcj->abci', iMM1a, full_data)

    #Move again axis to original positoin
    data_corrected=np.moveaxis(data_corrected,0,-1)
    data_corrected=np.moveaxis(data_corrected,0,-1)
    return data_corrected, MM1a
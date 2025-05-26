"""
This module contains the functions to correct for crosstalk
through the method described in Jaeggli et al. 2022. 
(https://doi.org/10.3847/1538-4357/ac6506)

Original functions:
https://github.com/sajaeggli/adhoc_xtalk/blob/main/Jaeggli_etal_2022ApJ_AdHoc_Xtalk.ipynb

v2: the Mueller matrix is fitted independently at each wavelength to take
into account spectral variations of the polarimetric response of the instrument
"""
import numpy as np
from matplotlib import colors, pyplot as plt
from astropy.io import fits
from glob import glob
from scipy.optimize import minimize
from scipy.ndimage import shift
from sklearn import linear_model
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
def fitfunc1(param, stokesin,wvl='all',method='jaeggli'):
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
    out = minimize_for_model1(iMM,stokesin,wvl=wvl,method=method)
    return(out)


def minimize_for_model1(iMM,bs,wvl='all',method='jaeggli'):
    """
    Function that computes the merit function that considers
    crosstalk from I to Q, U and V along the spectral profile
    for a given Mueller matrix
    Input:
        iMM: inverse of Mueller matrix of the diattenuator
        bs: measured Stokes I
        wvl: 'all' or selected wavelength sample (0,1,2,...).
            Wavelength to be employed to correct crosstalk.
        method: 'jaeggli', 'all_wvls' or 'corr'. Merit function
            ->jaeggli: metric defined in Eq. (16) of Jaeggli et al. (2022)
            ->all_wvls: modified Jaeggli's metric to compute the
                correlation between Stokes I at all wavelengths with
                Q, U and V at each wavelength
            ->corr: metric that computes the correlation of Stokes
                I with Q, U and V along the spectral profile as the
                 "sample correlation coefficient" defined in
                in https://en.wikipedia.org/wiki/Correlation
    Output:
        out: computed merit function            
    """
    # apply a mueller matrix (rotation) to a 3D stokes vector 
    # (npixel,wavelength,4)
    new_stokes = np.einsum('ij,abj->abi',iMM, np.squeeze(bs))
    Nwaves=new_stokes.shape[1] #Number of wavelenth samples

    # Minimization criteria
    if method=='jaeggli' or method=='corr':
        stI=new_stokes[:,:,0]
        stQ=new_stokes[:,:,1]
        stU=new_stokes[:,:,2]
        stV=new_stokes[:,:,3]
        if method=='jaeggli':
            kappa1=1
            kappa2=1
            kappa3=1
        elif method=='corr':
            #Subtract mean over spectral profile
            stI=stI-np.mean(stI,axis=1,keepdims=True)
            stQ=stQ-np.mean(stQ,axis=1,keepdims=True)
            stU=stU-np.mean(stU,axis=1,keepdims=True)
            stV=stV-np.mean(stV,axis=1,keepdims=True)

            #Compute normalization factors
            kappa1=np.sqrt(np.sum(stI**2,axis=1)*np.sum(stQ**2,axis=1)) 
            kappa2=np.sqrt(np.sum(stI**2,axis=1)*np.sum(stU**2,axis=1)) 
            kappa3=np.sqrt(np.sum(stI**2,axis=1)*np.sum(stV**2,axis=1))
    
        #Compute merit function
        if wvl=='all': #Correlate Q,U,V with I at each wavelength     
            out = np.abs(np.sum(stI*stQ,axis=1)/kappa1)\
                +np.abs(np.sum(stI*stU,axis=1)/kappa2)\
                +np.abs(np.sum(stI*stV,axis=1)/kappa3)
        else: # Correlate Q,U,V at selected wavelength with I at all wavelengths
            out = np.abs(stQ[:,wvl]*np.sum(stI,axis=1)/kappa1)\
                +np.abs(stU[:,wvl]*np.sum(stI,axis=1)/kappa2)\
                +np.abs(stV[:,wvl]*np.sum(stI,axis=1)/kappa3)    
    elif method=='all_wvls': #Correlate Q,U,V at all wavelengths with I at all wavelengths
        out1=0
        out2=0
        out3=0
        for i in range(Nwaves):
            for j in range(Nwaves):
                out1 += new_stokes[:,i,0]*new_stokes[:,j,3]
                out2 += new_stokes[:,i,0]*new_stokes[:,j,2] 
                out3 += new_stokes[:,i,0]*new_stokes[:,j,1]
        out=np.abs(out1)+np.abs(out2)+np.abs(out3)
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
                       last_wvl=None,plots=False,method='jaeggli'):
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
        method: 'jaeggli' (default), 'all_wvls'
            -> jaeggli: Jaeggli et al. (2022) modified merit function
            that minimizes the correlation between Stokes Q, U and V at
                indepdentently each wavelength with Stokes I at all wavelengths.
            -> all_wvls: modified merit function that minimizes altogether the
             correlation of Stokes Q, U and V at all wavelengths with Stokes I
             at all wavelengths 
    Output:
        datarest: 4D array with the corrected Stokes parameters
        MM1a: Diattenuation Mueller matrix that converts the
            "real" Stokes parameters into the "observed" ones.
        
    """
    x0=region[0]
    xf=region[1]
    y0=region[2]
    yf=region[3]
    #Reorder axis to convert into dimensions: [x,y,wavelength,stokes]
    data=np.moveaxis(data,0,-1)
    data=np.moveaxis(data,0,-1)

  
    if norm is True:
        #Normalization of data
        norm_factor=np.median(data[:,:,0,0])
        data=data/norm_factor


    ## Crop data to avoid edge effects arising from alignment/rotation
    full_data=data.copy()
    data=data[x0:xf,y0:yf,:,:]


    #Last wavelength to be considered in the minimization
    data=data[:,:,:last_wvl,:]
    Nwaves=data.shape[2]

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
        axs[0,0].set_title('Stokes I (cont)')
        axs[0,1].imshow(data[:,:,0,1],cmap='gray')
        axs[0,1].set_title('Stokes Q (cont)')
        axs[1,0].imshow(data[:,:,0,2],cmap='gray')
        axs[1,0].set_title('Stokes U (cont)')
        axs[1,1].imshow(data[:,:,Nwaves//2,3],cmap='gray')
        axs[1,1].set_title('Stokes V (core)')
        for i in range(2):
            for j in range(2):
                axs[i,j].contour(pmap, [pthresh], colors='green', 
                                linewidths=0.75)
        plt.show()
        plt.close() 

    #Fit the Mueller matrix at each wavelength
    data_corrected=full_data.copy() 
    MM1a=np.zeros((Nwaves,4,4))   
    for wvli in range(Nwaves):
        #Wrap function to use it in minimize with positional and keyword arguments
        fun=lambda x: fitfunc1(x, weak_region, wvl=wvli, method=method)
            
        #Minimize merit function
        result = minimize(fun, initial_guess)

        # Apply correction for I<->QUV cross-talk at each wavelength
        MM1a[wvli,:,:] = polmodel1(result.x[0],result.x[1], result.x[2])
        iMM1a = np.linalg.inv(MM1a[wvli,:,:])
        data_inverted =  np.einsum('ij,abcj->abci', iMM1a, full_data)
        data_corrected[:,:,wvli,:] = data_inverted[:,:,wvli,:]

    #Move again axis to original positoin
    data_corrected=np.moveaxis(data_corrected,0,-1)
    data_corrected=np.moveaxis(data_corrected,0,-1)
    return data_corrected,MM1a


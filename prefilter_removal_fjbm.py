"""
This module fits the prefilter using the blueshift
map and the gain map obtained with the program "compute_master_flat_field.py".

We correct from blueshift the spectral profile of the flat while
assuming that the spectral positions of the FTS are correct and do not need
to be shifted

#
"""   
import sys
import os
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
sys.path.append('./functions')
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import interp1d,RectBivariateSpline
#from scipy.signal import convolve
from numpy import convolve

"""
Input parameters
"""
#Etalon parameters (only to compute the effective reflectivity)
n=2.5 #As fitted by Pablo
h=281e-6 #As fitted by Pablo


def gaussian(x,x0,sigma):
    """
    Gaussian function defined so FWHM=2.355*sigma
    Input:
        x: independent variable
        x0: center of the Gaussian
        sigma: standard deviation
    """
    return np.exp(-0.5*((x-x0)/sigma)**2)

def prefilter_model(om,wvlv):
    """
    Function to create a prefilter model depending on the observation mode
    based on the fitted parameters for the prefilter.
    """
    if om=='2.02':
        pref_wvl=525.04e-9
        pref_sigma=0.055e-9
        wvl0=525.0217e-9-0.125e-12 #Central wavelength of the line (according to FTS)
    elif om=='2.06':
        pref_wvl=525.07e-9
        pref_sigma=0.04e-9
        wvl0=525.0653e-9   
    elif om=='1':
        pref_wvl=517.26e-9
        pref_sigma=0.05e-9
        wvl0=517.27e-9
    pref_model=gaussian(wvlv+wvl0,pref_wvl,pref_sigma)
    return pref_model


def prefilter_fitting(cam_ave,om,wvlv):
    """
    Function to fit the prefilter parameters over the mean profile of the flats
    """   
    #Central wavelength of the observed line and optimization parameters
    meth='Nelder-Mead' #Minimization method

    if om=='2.02':
        wvl0=525.0217e-9-0.125e-12 #Central wavelength of the line (according to FTS)
        init_guess=[3e-12,525.05e-9,0.06e-9] #Initial guess [etalon width (m), pref.  wvl (m), pref. width (m)]
    elif om=='1':
        wvl0=517.27e-9
        init_guess=[3.5e-12,517.25e-9,0.06e-9] #Initial guess [etalon width (m), pref.  wvl (m), pref. width (m)]
    elif om=='2.06':
        wvl0=525.0653e-9
        init_guess=[3.8e-12,525.05e-9,0.06e-9]
    bounds=[(1e-12,10e-12), #Transmission width
            (wvl0-1e-9,wvl0+1e-9), #Central wavelength
            (0.01e-9,1e-9)] #Prefilter width


    """
    Interpolate FTS
    """
    #Import FTS
    file = os.path.join(__location__, 'fts.npz')
    data = np.load(file)
    fts0 = data['fts']/10000 #Normalize FTS spectrum to 1
    fts_w = data['fts_w']*1e-10

    #Interpolate FTS over wavelength
    fts_interp=interp1d(fts_w,fts0,kind='linear')

    #Convolve FTS with different Gaussians and interpolate in 2D
    sigma_vector=np.arange(1e-12,11e-12,1e-12) #Width of the Gaussian in meters
    wvlv_interp=np.linspace(-0.1e-9,0.1e-9,10001) #Interval of wvls for interpolation centered about 0

   

    #Loop over different sigmas to create 2D interpolation
    i=-1
    fts_matrix=np.zeros((len(wvlv_interp),len(sigma_vector)))
    for sigma in sigma_vector:
        i+=1
        norm=np.sum(gaussian(wvlv_interp,0,sigma))
        fts_matrix[:,i]=convolve(fts_interp(wvl0+wvlv_interp),
                            gaussian(wvlv_interp,0,sigma),mode='same')/norm
        #Shift to center minimum
        min_wvl_conv=minimize(lambda x: fts_matrix[:,i][np.argmin(np.abs(wvlv_interp - x))],
                            x0=0,method='Nelder-Mead',bounds=[(-1e-12,1e-12)]).x[0]
    
        fts_matrix[:,i]=np.interp(wvlv_interp + min_wvl_conv,wvlv_interp,fts_matrix[:,i])
    fts_interp2D=RectBivariateSpline(wvl0+wvlv_interp,sigma_vector,fts_matrix,kx=1,ky=1)


    """
    Shift FTS to match the observed central wavelength
    """
    #Compute minimum of the average profile over x, y and the 4 modulations
    cam_ave_interp=interp1d(wvlv,cam_ave,kind='quadratic',bounds_error=False,
                            fill_value='extrapolate')

    wvlv2=np.linspace(wvlv[0],wvlv[-1],5000)
    min_index=np.argmin(cam_ave_interp(wvlv2))
    wvl_min= wvlv2[min_index]


    #Normalization factor to match the intensity of the profiles with FTS
    if om=='2.02':
        norm=1*cam_ave[-1]/fts_interp(wvl0+wvlv[-1])
        prof_scaled=cam_ave/norm
    elif om=='2.06': 
        norm=1.085*cam_ave[-1]/fts_interp(wvl0+wvlv[-1])
        prof_scaled=cam_ave/norm 
    elif om=='1':
        norm=1.15*cam_ave[0]/fts_interp(wvl0+wvlv[0])
        prof_scaled=cam_ave/norm


    """
    Fit the transmission and prefilter parameters over the region with 
    minimum mismatch
    """
    #Shift needed for the profile to coincide with FTS minimum
    delta_wvl=wvl0 - wvl_min

        
    #Convolve FTS with Gaussian to match the resolution
    def merit_indiv1(params):
        sigma=params[0] #Sigma of transmission profile
        wvl_pf=params[1] #Central wavelength of prefilter
        pf_width=params[2] #Sigma of prefilter

        #Shifted wavelength axis
        wvlv_shifted=wvlv + delta_wvl

        #Interpolated mean profile at shifted wavelengths
        fts_shifted=fts_interp2D(wvlv_shifted, sigma).flatten()

        #Compute difference with individual profile avoiding continuum wavelength
        diff=prof_scaled[:] - fts_shifted[:]*gaussian(wvlv_shifted,wvl_pf,pf_width)  
        return np.sum(diff**2)


    #Minimize merit function
    minim=minimize(merit_indiv1,init_guess,method=meth,bounds=bounds)#,options=opt)
    sigma_opt=minim.x[0] #Optimized etalon width
    wvl_pf_opt=minim.x[1]+wvl_min #Correct prefilter wvl for shift
    pf_width_opt=minim.x[2]


    print('Fitted prefilter parameters for om='+om+':')
    print('Etalon sigma (pm): ',np.round(sigma_opt*1e12,2))
    print('Prefilter central wavelength (nm): ',np.round(wvl_pf_opt*1e9,3))
    print('Prefilter sigma (nm): ',np.round(pf_width_opt*1e9,3))

    #Compute prefilter
    prefilter=gaussian(wvlv+wvl0,wvl_pf_opt,pf_width_opt)
    return prefilter


def etalon_width(R,wvl0,n,h):
    #Compute etalon width given reflectivity R
    finesse=np.pi*np.sqrt(R)/(1-R)
    return wvl0**2/(2*n*h*finesse)

def effective_reflectivity(fwhm_opt,wvl0,n,h):
    #Find effective reflectivity that reproduces the observed FWHM
    R_eff=minimize(lambda R: abs(etalon_width(R,wvl0,n,h)-fwhm_opt)**2,0.92,
                   method='Nelder-Mead',bounds=[(0.7,0.95)])
    return R_eff.x[0]


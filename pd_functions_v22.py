"""
Same as pd_functions_v16, except for the fact that now admits
a phase diversity different from a defocus. In particular,
the parameter a_d can be now an array. If so, it must contain
the Zernike coefficients of the phase diversity in order.

Changes with respect to pd_functions_v16
-'pupil' and 'wavefront' now differentiate between a_d being a scalar or an
 array. If a_d is an array, it is interpreted as the Zernike coefficients
 introduced by the phase diversity image.
-Normalization is corrrect in read_image, so the focused and defocused
images are normalized independently, each one by their own mean.
-low_f is increased up to 0.2 and a second low filtering is applied
to the noise filter in order to elliminate isolated peaks

Changes with respect to pd_functions_v17
-prepare_PD now uses the Sicairos method (dftreg) to
align the FFTs with subpixel accuracy before feeding
the PD algorithm (modification made on 2023/3/10)
-we compute the optimized merit function when calling to "object_estimate"
(modificatino made on 2023/3/16)
-noise_filter is calculated in the 4 corners of the FFT image
-gamma is calculated with non-apodized images, like in v1.6
-A correction in filter_sch is applied. The denominator was calculated
using gamma1(=1) instead of the calculated gamma.
-Tip/tilt correction is applied ONLY to defocused images in case Jmin<4

Changes with respect to pd_functions_v20
-Tip/tilt correction is applied individually to each OTF. Note that in
pd_functions_v20 tip/tilt was the same for all OTFs. Changes affect OTF.py, which
interprets the first 2K+1 terms as a global offset and 2K tip/tilt coefficients,
each one corresponding to each image (0 for the focused image)
-np.int substituted by "int" as it was deprecated in NUmpy 1.20
-We introduce two regularization factors in Qfactor (reg1*(nu/nuc)^reg2) for
reconstruction purposes. Default values: reg1=0.02,reg2=1

Changes with respect to pd_functions_v21
- We introduce a new merit function that takes into account a jitter term
to be adjusted between two images: one with and another one without jitter.
- Default parameters (N, R, nuc, etc.) are no longer set as default, but
introduced by the user or computed within a function 
"""
#from skimage.morphology import square, erosion
from matplotlib import pyplot as plt
from PIL import Image as im
from photutils import CircularAperture #conda install -c astropy photutils
from scipy.ndimage import median_filter
from scipy import ndimage
from scipy.fftpack import fftshift, ifftshift, fft2, ifft2
from scipy.optimize import minimize, minimize_scalar
from astropy.convolution import convolve
import numpy as np
import zernike as zk
import math_func2 as mf
import sys
import os
import time
import scipy
import shift_func as sf
from astropy.io import fits
try:
    import pyfftw
    flag_pyfftw=1
except:
    flag_pyfftw=0


#TuMag parameters



def tumag_params():
    wvl=525.02e-9 #517.3e-9 #Wavelength [m]
    fnum=60 # f-number
    Delta_x=11e-6 ##Size of the pixel [m]
    return wvl,fnum,Delta_x

def compute_nuc(N,wvl,fnum,Delta_x):
    """
    This function computes the critical frequency (in pixel units)
    for an image of size NxX. It returns also the pupil radius
    (in pixel units)
    """
    nuc=1/(wvl*fnum) #Critical frequency (m^(-1))
    inc_nu=1/(N*Delta_x)
    R=(1/2)*nuc/inc_nu #Pupil radius [pixel]
    nuc=2*R #critical frequency (pixels)
    nuc-=1 #To avoid boundary zeros when calculating the OTF (Julian usa 1
    return nuc,R

def sampling(N):
    """
    This function creates a grid of points with NxN dimensions for calling the
    Zernike polinomials.
    Output:
        X,Y: X and Y meshgrid of the detector
    """

    if N%2 != 0: #If N is odd (impar)
        print('Error in function "sampling": Number of pixels must be even!')
        x=np.linspace(-round(N/2),round(N/2),N)
        #sys.exit()
    else: #If N is even (par)
        x=np.linspace(-N/2,N/2,N)
        #x=np.linspace(-N/2,N/2,N+1)
        #x=x[:-1]

    y=np.copy(x)
    X,Y=np.meshgrid(x,y)

    return X,Y

def sampling2(N,R):
    """
    It returns RHO and THETA with the sampling of the pupil carried out
    in fatima.py. RHO is already normalized to the unity.
    """
    r=1
    x=np.linspace(-r,r,2*int(R))
    y=np.copy(x)
    [X,Y]=np.meshgrid(x,y)
    RHO = np.sqrt(X**2+Y**2)
    THETA = np.arctan2(Y, X)
    RHO0=np.zeros((N,N))
    THETA0=np.copy(RHO0)
    RHO0[N//2-int(R):N//2+int(R),N//2-int(R):N//2+int(R)]=RHO
    THETA0[N//2-int(R):N//2+int(R),N//2-int(R):N//2+int(R)]=THETA
    return RHO0,THETA0

def cir_aperture(R,N,ct=None):

    #Offset to be subtracted to center the circle
    if ct==None:
        if N%2 != 0: #If N is odd (impar)
            #print('Number of pixels should be an even integer!\n')
            ct = 1/2
        else: #If N is even (par)
            ct = -1/2#1/2

    N = int(N)
    R = int(R)
    A = CircularAperture((N/2-ct,N/2-ct),r=R) #Circular mask (1s in and 0s out)
    A = A.to_mask(method='exact') #Mask with exact value in edge pixels
    #A = CircularAperture((N/2-ct,N/2-ct),r=R-1)
    #A = A.to_mask(method='center')

    A = A.to_image(shape=(N,N)) #Conversion from mask to image
    return A



def pmask(nuc,N):
    """
    Mask equal to 1 where cir_aperture is 0
    """
    pmask=np.where(cir_aperture(nuc,N,ct=0)==0,1,0)
    return pmask

def aperture(N,R,cobs=0,spider=0):
    """
    This function calculates a simple aperture function that is 1 within
    a circle of radius R, takes and intermediate value between 0
    and 1 in the edge and 0 otherwise. The values in the edges are calculated
    according to the percentage of area corresponding to the intersection of the
    physical aperture and the edge pixels.
    http://photutils.readthedocs.io/en/stable/aperture.html
    Input:
        N: 1D size of the detector
        R: radius (in pixel units) of the pupil in the detector
        cobs: central obscuration (as percentage of radius)
        spider: width of spider arms (in pixels)
    Output:
        A: 2D array with 0s and 1s
    """
    A=cir_aperture(R,N)

    #If central obscuration:
    if (cobs != 0):
        if N%2 != 0:
            #print('Number of pixels should be an even integer!\n')
            ct = 0.
        else:
            ct = 1/2
        B=CircularAperture((N/2-ct,N/2-ct),r = R*cobs/100.)
        B=B.to_mask(method='exact') #Mask with exact value in edge pixels
        B=B.to_image(shape=(N,N)) #Conversion from mask to image
        A=A-B
        A = np.where(A<=0,0,A)

    #If spider:
    if (spider != 0):
        C=array((N*10,N*10))
        M = N*10 + 10
        S = R*spider/10.
        C[int(M/2-S/2):int(M/2+S/2),:] = 1.0
        C[:,int(M/2-S/2):int(M/2+S/2)] = 1.0
        nC = imresize(C, (N, N))#, interp="bilinear")
        nC = nC / np.max(nC)
        nC=np.where(nC<=0,0,nC)
        A = A - nC
        A = np.where(A<=0,0,A)
    return A

def aperture_easy(X,Y):
    """
    This function calculates a simple aperture function that is 1 within
    a circle of radius R and 0 otherwise.
    Input:
        X,Y: meshgrid with the coordinates of the detector ('sampling.py')
        R: radius (in pixel units) of the pupil in the detector
    Output:
        A: 2D array with 0s and 1s
    """
    A=X**2+Y**2<=R**2
    return A

def pupil(a,a_d,RHO,THETA,A):
    """
    This function calculates the generalized pupil function of a circular
    aperture for an incident wavefront with an aberration given by the
    Zernike coefficients following Noll's order.
    The aperture function is given by A ('aperture.py').
    Input:
        a: 1D array with Zernike coefficients (ordered following Noll's order
           The first element,a[0], is the Piston, a[1] and a[2] are Tilts
           in the X and Y direction ...), except if tiptilt is True
        a_d: Defocus Zernike coefficient for defocused image (float or int)
            or array with ALL the aberrations introduced by the defocused image
        RHO,THETA: meshgrid with polar coordinates (rho normalised to 1)
        A: aperture 2D array ('aperture')
    Output:
        p: 2D complex array describing the pupil function
    """
    try:
        Jmax=len(a) #Max. number of Zernikes to be employed
    except TypeError:
        print('Error in "pupil" function: "a" must be a Numpy array or a list')
        sys.exit()

    #Phase induced during PD
    if isinstance(a_d, (float,int)): #If float, then induce a defocus a_d
        phi=a_d*zk.zernike(0,2,RHO,THETA)#diversity phase
    #If vector, a_d contain all the aberrations of the defocused image
    else:
        print('WARNING: a_d interpreted as a vector with the aberrations of the\
        PD plate in "pupil.py"')
        #If array, then get zernike coefficients from a_d
        Jmax_div=len(a_d)
        jj=0
        phi=0*RHO #To have an array filled with zeroes
        while jj<Jmax_div:
            jj+=1
            phi+=a_d[jj-1]*zk.zernikej_Noll(jj,RHO,THETA)

    #Wavefront error produced by Zernike ('a') coefficients
    if a_d==0: #If focused image, then do not consider offset and tip/tilt
        a[:3]=0

    j=0
    while j<Jmax:
        j+=1
        phi+=a[j-1]*zk.zernikej_Noll(j,RHO,THETA)
    phi=A*phi
    p=A*np.exp(1j*phi)
    return p

def wavefront(a,a_d,RHO,THETA,A,R,N):
    """
    This function returns the wavefront map given a set of Zernike
    coefficients. It works similarly than pupil.py
    """
    try:
        Jmax=len(a) #Max. number of Zernikes to be employed
    except TypeError:
        print('Error in "pupil" function: "a" must be a Numpy array or a list')
        sys.exit()
    #Phase diversity
    if isinstance(a_d, (float,int)): #If float, then induce a defocus a_d
        phi=a_d*zk.zernike(0,2,RHO,THETA)#diversity phase
    else: #If array, then get zernike coefficients from a_d
        Jmax_div=len(a_d)
        jj=0
        phi=0*RHO #To have an array filled with zeroes
        while jj<Jmax_div:
            jj+=1
            phi+=a_d[jj-1]*zk.zernikej_Noll(jj,RHO,THETA)
    j=0
    while j<Jmax:
        j+=1
        phi+=a[j-1]*zk.zernikej_Noll(j,RHO,THETA)
    phi=A*phi #Calcular phi en todo X,Y para luego hacer esto
                              #puede no ser lo más óptimo
    #Subframing to cover only the pupil
    xc=int(N/2)
    x0=int(xc-R+1)
    xf=int(xc+R+1)

    phi=phi[x0:xf,x0:xf]

    return phi

def PSF(a,a_d,RHO,THETA,ap):
    """
    This function calculates the PSF of a circular aperture for an incident
    wavefront with a wavefront aberration given by the Zernike coefficients.
    The aperture function is defined in 'aperture.py'
    Input:
        a,a_d,RHO,THETA,A: defined in 'pupil' function
    Output:
        psf: 2D real array representing the PSF
    """
    #Generalized complex pupil function
    p=pupil(a,a_d,RHO,THETA,ap)
    #PSF
    psf=np.abs(mf.fourier2(p))**2 #Non normalized PSF (normalization made in OTF)
    return psf


def OTF(a,a_d,RHO,THETA,ap,norm=None,K=2,tiptilt=False):
    """
    This function calculates the OTFs of a circular aperture for  incident
    wavefronts with aberrations given by a set of Zernike coefficients.
    The OTF is calculated as the autocorrelation of the pupil.
    Input:
        a,a_d,RHO,THETA,A: defined in 'pupil.py'
            a_d can be an array with shape K, where K is the number
            of images with different defocuses. In such a case, a_d
            contains the defocus coefficient of each image and the program
            returns an multidimensional array of OTFs
        norm:{None,True}, optional. 'True' for normalization purpose. Default
         is None
        tiptilt: {True,False). If True, the first 2*(K-1) Zernike terms
        correspond to tip/tilt terms.
    Output:
        otf: 2D complex array representing the OTF if len(a_d)!=K or
            a 3D complex array whose 3rd dimenstion indicates the  OTF of
            of each defocused image.
        norm: normalization factor of the OTF

    """
    #If a_d is a number or a vector containing the aberrations of the PD plate
    if isinstance(a_d, (float,int)) or len(a_d)>K:
        #Pupil
        p=pupil(a,a_d,RHO,THETA,ap)
        #OTF
        norma_otf,otf=mf.corr(p,p,norma=True)
        if norm==True:
            norma=norma_otf
            #norma=np.max(np.abs(otf)[:])
            otf=otf/norma #Normalization of the OTF
        else:
            norma=1
        otf=otf[...,np.newaxis]#To create a 3rd dummy axis    
    #If a_d is an array containing K diversities
    elif len(a_d)==K:
        otf=np.zeros((RHO.shape[0],RHO.shape[0],K),dtype='complex128')
        norma=np.zeros(K)

        for i in range(K):
            #Select tiptilt terms for each image along the series
            if tiptilt is True:
                #Offset and tiptilt terms corresponding to each position
                a1=select_tiptilt(a,i,K)
            else:
                a1=a
            #Pupil computation
            p=pupil(a1,a_d[i],RHO,THETA,ap)

            #OTF
            norma_otf,otf[:,:,i]=mf.corr(p,p,norma=True)
            if norm==True:
                norma[i]=norma_otf
                #norma=np.max(np.abs(otf)[:])
                otf[:,:,i]=otf[:,:,i]/norma[i] #Normalization of the OTF
            else:
                norma[i]=1
    return otf,norma

def OTF_jitt(sigma,nuc,N,norm=None):
    """
    This function computes the OTFs of jitter corresponding to a pair of images:
    one free of jitter and another one affected by it. It is based on the
    model of Fiete 2014 (Eq. 19)
    Input:
        sigma: vector with the rms of jitter in x,y and its correlation
            (sigma_x,sigma_y,sigma_xy)
            
        norm:{None,True}, optional. 'True' for normalization purpose. Default
         is None
        
    Output:
        otf: 2D complex array representing the OTFs of the image free of
        jitter (all ones), and of the image affected by it.
        norm: normalization factor of the OTF

    """
    #Meshgrid with frequency coordinates
    delta_theta=0.0378 #Plate scale in arcseconds (arcsec/pixel)
    nu_lim=1/(2*delta_theta) 
    nu_vec=np.linspace(-nu_lim,nu_lim,int(N)+1) #Odd number of samples to center the stray-light OTF
    NU,ETA=np.meshgrid(nu_vec,nu_vec)
    NU=NU[:-1,:-1] #To have again an even number of samples
    ETA=ETA[:-1,:-1]

    #OTF_jitter model for the jitter-free and the jitter-affected images
    otf_jitt=np.zeros((N,N,2))
    otf_jitt[:,:,0]=np.ones((N,N))
    otf_jitt[:,:,1]=np.exp(-2*np.pi**2*(NU**2*sigma[0]**2+ETA**2*sigma[1]**2+NU*ETA*sigma[2]**2))
    
    if norm==True:
        norma=np.max(np.abs(otf_jitt)[:])
        otf=otf/norma #Normalization of the OTF
    else:
        norma=1
    
    return otf_jitt,norma

def select_tiptilt(a,i,K):
    """
    Returns Zernike vector (a) with selected tip tilt terms for each position
    (i) along the series of K images
    """
    tiltx=float(a[(2*i+1)])
    tilty=float(a[(2*i+2)])
    firsta=np.zeros((3,1))
    firsta[0]=0
    firsta[1]=tiltx
    firsta[2]=tilty
    a1=np.concatenate((firsta,a[(2*K+1):])) #0 is for the offset term
    return a1


def OTF_circular(nuc,N):
    """
    This function returns the analytical expression of the OTF of a circular
    aperture with no aberrations.
    Output:
        otf: 2D complex array representing the OTF of a circular aperture
    """
    nu=np.linspace(-N/2,N/2,N)
    nu_norm=nu/nuc
    mask=np.abs(nu_norm**2)<=1
    otfc_out=np.zeros(N)
    otfc_in=(2/np.pi)*(np.arccos(np.absolute(nu_norm*mask))-\
    np.absolute(nu_norm*mask)*np.sqrt(1-(nu_norm*mask)**2))
    otfc=np.where(mask,otfc_in,otfc_out)
    return otfc

def convPSF(I,a,a_d,RHO,THETA,ap,norm=None):
    """
    This function calculates the convolution of an image with the PSF
    of IMaX+. 'I' must be a Numpy array of type 'float64'
    Input:
        I: (real) 2D numpy array (image)
        a,a_d,RHO,THETA,A: defined in 'pupil.py'
    Output:
        d: real 2D numpy array representing the convolution of I and the system
    """
    if I.dtype != 'float64':
        print(I.dtype,'Error: I must be a Numpy array of type float64')
        sys.exit()
    #Fourier transform of the image
    O=fft2(I)
    O=ifftshift(O)

    #Aberrated image
    otf,norma=OTF(a,a_d,RHO,THETA,ap,norm=norm)
    D=O*otf

    d=ifftshift(D)
    d=ifft2(d).real
    return d

def Zkfactor(k,Ok,Hk,Q,gamma):
    """
    Factor Zk defined in Paxman 1992 (eq. C6) needed to calculate
    the derivative of the merit function (eqs. C1-C5)
    Input:
        k:index 'k' of Zk (k<= #of different defocuses)
        Ok: 3D array with the FFTs of each of the PD images with different defoc.
        Hk: 3D array with the OTFs of the system for each defocus
    Output:
        Zk: 2D array with the Zk factor
    """
    Qinver=Qinv(Hk,gamma)
    sum1=Qinver**2*np.sum(Ok*np.conj(Hk),axis=2)*np.conj(Ok[:,:,k])
    sum2=np.abs(np.sum(Ok*np.conj(Hk),axis=2))**2*gamma[k]*np.conj(Hk[:,:,k])
    Zk=Q**4*(sum1-sum2)
    return Zk


def Qfactor(Hk,gamma,nuc,N,reg1=0.02,reg2=1):
    """
    Q factor defined in Lofdahl and Scharmer 1994 for construction of the merit
    function
    Input:
        Hk: 3D array with the OTFs of the system for each defocus
    Output:
        Q: 2D complex numpy array representing the Q factor
    """
    np.seterr(divide='ignore')
    #For the effect of the regularization to increase with the frequency in the way
    #described in Martinez Pillet (2011), Sect. 9.2)
    if reg1>0:
        nx=Hk.shape[0]
        ny=Hk.shape[1]
        x=np.linspace(-int(nx/2),int(nx/2),nx)
        y=np.linspace(-int(ny/2),int(ny/2),ny)
        [X,Y]=np.meshgrid(x,y)
        Q=1/np.sqrt(np.sum(gamma*np.abs(Hk)**2,axis=2)+reg1*(np.sqrt(X**2+Y**2)/nuc)**reg2)
    else:
        Q=1/np.sqrt(np.sum(gamma*np.abs(Hk)**2,axis=2)+reg1)
    
    Q=np.nan_to_num(Q, nan=0, posinf=0, neginf=0)
    Q=Q*cir_aperture(nuc,N,ct=0)
    return Q


def Qinv(Hk,gamma,nuc,N):
    """
    Inverse of Q. Defined for test purposes
    Input:
        Hk: 3D array with the OTFs of the system for each defocus
    Output:
        Qinv:
    """
    Qinv=np.sqrt(np.sum(gamma*np.abs(Hk)**2,axis=2))
    Qinv=Qinv*cir_aperture(nuc,N,ct=0)
    return Qinv

def Ffactor(Q,Ok,Hk,gamma):
    """
    F factor defined by Lofdahl and Scharmer 1994 (Eq. 5) in the general
    form of Paxman 1992 (Eq. 19). Gamma is added, too.
    Input:
        Q: Q returned in Qfactor
        Ok: 3D array with the FFTs of each of the PD images with different defoc.
        Hk: 3D array with the OTFs of the system for each defocus
    Output:
        Zk: 2D array with the Zk factor
    """

    return Q**2*np.sum(gamma*Ok*np.conj(Hk),axis=2)

def smooth_filt(array,size=3):
    """
    Function that smoothes a 2D array using ndimage.uniform_filter
    """
    #kernel = np.ones((size,size))
    #return convolve(array, kernel)
    return ndimage.uniform_filter(array, size=size)


def noise_power(Of,nuc,N,filterfactor=1.5):
    """
    Average level of noise power of the image. Based on Bonet's program
    'noise_level.pro'. Noise is calculated on 4 quadrants of the
    FFT of the image, beyond the critical frequency, to elliminate
    the dreadful cross (horizontal and vertical spurious signal in the FFT)
    that appears because of the finite Nyquist frecuency of the detector.
    """
    #Circular obscuration mask to calculate the noise beyond the critical freq.
    cir_obs=pmask(nuc,N)


    #Calculation of noise
    #power=np.sum(np.abs(Of)**2*cir_obs)/np.sum(cir_obs)

    #1st quadrant
    x2=int(np.floor(N/2-nuc*np.sqrt(2)/2))
    power=np.sum((np.abs(Of)**2*cir_obs)[0:x2,0:x2])/np.sum(cir_obs[0:x2,0:x2])
    #2nd quadrant
    x3=N-x2
    power+=np.sum((np.abs(Of)**2*cir_obs)[x3:N,0:x2])/np.sum(cir_obs[x3:N,0:x2])
    #3rd quadrant
    power+=np.sum((np.abs(Of)**2*cir_obs)[0:x2,x3:N])/np.sum(cir_obs[0:x2,x3:N])
     #4th quadrant
    power+=np.sum((np.abs(Of)**2*cir_obs)[x3:N,x3:N])/np.sum(cir_obs[x3:N,x3:N])

    #To obtain a more conservative filter in filter_sch
    power=filterfactor*power
    return power

def filter_sch(Q,Ok,Hk,gamma,nuc,N,low_f=0.2):
    """
    Filter of the Fourier transforms of the focused and defocused
    object (Eqs. 18-19 of Lofdahl & Scharmer 1994). Based on
    filter_sch.pro (Bonet's IDL programs)
    Input:
        Q: Q returned in Qfactor
        Ok: 3D array with the FFTs of each of the PD images with different defoc.
        Hk: 3D array with the OTFs of the system for each defocus
    Output:
        filter:2D array (float64) with filter
    """
    denom=np.abs(Ffactor(Q,Ok,Hk,gamma))**2\
    *Qinv(Hk,gamma,nuc,N)**2

    filter=noise_power(Ok[:,:,0],nuc,N)/smooth_filt(denom)
    filter=(filter+np.flip(filter))/2
    filter=1-filter
    filter=np.where(filter<low_f,0,filter)
    filter=np.where(filter>1,1,filter)
    #filter=erosion(filter, square(3)) #To remove isolated peaks (also decreases
                                       #the area of the central region)

    filter=median_filter(filter,size=9)
    filter=smooth_filt(filter)
    filter=filter*cir_aperture(nuc,N,ct=0)
    filter=np.nan_to_num(filter, nan=0, posinf=0, neginf=0)

    #Apply low_filter again to elliminate isolated peaks
    filter=np.where(filter<low_f,0,filter)
    return filter


def meritE(Ok,Hk,Q):
    """
    Merit function 'E' for phase diversity optimization
    Ref: Lofdahl and Scharmer 1994, Eq.9
    Input:
        Q: 2D array returned by Qfactor
        Ok: 3D array with the FFTs of each of the PD images with different defoc.
        Hk: 3D array with the OTFs of the system for each defocus
    Ouput:
        E: 2D complex array with E factor for construction of the merit function
            or 3D complex array with the error metrics of each possible
            pair of the K images
        function
    """
    #Lofdahl error metric
    #E=Q*(Ok[:,:,1]*Hk[:,:,0]-Ok[:,:,0]*Hk[:,:,1])

    #General error metrix from Paxman merit function
    #E=Q*np.sum(Ok*np.conj(Hk),axis=2)

    #Alternative general error metric
    #E=Q*sumamerit(Ok,Hk)


    #Array with invidual error metrics
    E=arraymerit(Ok,Hk)
    if E.ndim==3:
        for i in range(E.shape[2]):
            E[:,:,i]=Q*E[:,:,i]
    elif E.ndim==2:
        E=Q*E
    return E

def derivmerit(Ok,Hk,derHk,Q,gamma,cut=None):
    """
    Derivative of the merit function of Paxman 1992 (eq. C5)
    Input:
        Q: 2D array returned by Qfactor
        Ok: 3D array with the FFTs of each of the PD images with different defoc.
        Hk: 3D array with the OTFs of the system for each defocus
        derHk: 3D array with the derivatives of the OTF for each defocus
    """


    sumZkSk=np.zeros((Ok.shape[0],Ok.shape[1]),dtype='complex128') #Initialization of the sum
    for k in range(Ok.shape[2]):
        sumZkSk+=Zkfactor(k,Ok,Hk,Q,gamma)*derHk[:,:,k]#1st term of eq. C5


    #Image domain
    sumZkSk=ifftshift(sumZkSk)
    N=sumZkSk.shape[0]
    sumzksk=N*N*ifft2(sumZkSk)
    sumzksk-=np.mean(sumzksk)

    #Parseval theorem
    n=sumzksk[cut:-cut,cut:-cut].shape[0]
    deriv=sumzksk[cut:-cut,cut:-cut]
    #deriv+=np.conj(deriv) #Second term of eq. C5

    deriv=2*np.real(deriv)#1st + 2nd term = 2*real(1st term)
    deriv=np.sum(deriv)/n**2

    return deriv

def sumamerit(A,B):
    """
    This function returns the double summatory from j=1 to K-1
    and from k=j+1 to K defined by Paxman (Eq. B10)to calculate the merit
    function.
    Input:
        A,B: 3D arrays of shape (N,N,K) containing Ok and Hk/derHk
    """
    K=A.shape[2] #Number of diversities (including focused image)
    suma=0 #Initialization of E
    for j in range(K-1):
        for k in range(j+1,K):
            suma+=A[:,:,j]*B[:,:,k]-A[:,:,k]*B[:,:,j]
    return suma

def arraymerit(A,B):
    """
    This function returns calculates the individual error
    functions of Eq. 20 in Paxman and it stores them in
    an array of shape (N,N,nk), where nk=(K-1)*K/2.
    Input:
        A,B: 3D arrays of shape (N,N,K) containing Ok and Hk or derHk
    """
    K=A.shape[2] #Number of diversities (including focused image)
    nk=int((K-1)*K/2)
    arraymerit=np.zeros((A.shape[0],A.shape[1],nk),dtype='complex128') #Initialization of E
    i=-1
    for j in range(K-1):
        for k in range(j+1,K):
            i+=1
            arraymerit[:,:,i]=A[:,:,j]*B[:,:,k]-A[:,:,k]*B[:,:,j]
    return arraymerit

def pair_index(i,K):
    m=-1
    for j in range(K-1):
        for k in range(j+1,K):
            m+=1
            if i==m and j==0:
                index=m
    return index

def meritL(E):
    """
    Merit function "L" defined by Löfdahl & Schamer 1994 (Eq. 8).
    """
    L=np.sum(np.abs(E)**2)#-np.sum(np.sum(np.abs(Ok)**2,axis=2))
    return L

def merite(E):
    """
    Inverse transform of 'E' ('meritE' function).
    If E is a 3D array, it returns also a 3D array with the IFFT of each metric.
    Input:
        E: 2D complex array with E factor for construction of the
            classical (Löfdahl & Scharmer) merit function
            or 3D complex array with the error metrics for each possible
            combination of pairs of the K images
    Output:
        e: IFFT of E (2D array or 3D depending of the shape of E)
    """
    if E.ndim==3: #Array with individual merit functions of K diverities
        e=np.zeros(E.shape)
        for i in range(E.shape[2]):
            e0=ifftshift(E[:,:,i])
            e0=ifft2(e0)
            n=E.shape[0]
            e0=n*n*e0.real
            e[:,:,i]=e0-np.mean(e0)#We substract the mean value (Bonet e_func.pro (CORE2007))
    elif E.ndim==2:
        e=ifftshift(E)
        e=ifft2(e)
        n=E.shape[0]
        e=n*n*e.real
        e-=np.mean(e)#We substract the mean value (Bonet e_func.pro (CORE2007))
    return e

def meritl(e,cut=None):
    """
    Merit function "L" defined by Lofdahl and Scharmer 1994.
    If 'e' is a 3D array, the merit function is generalized to work for
    K images with different defocuses
    """
    cutx=cut
    cuty=cut
    L=0
    if e.ndim==3: #If dim=3: sum individual merit functions
        #print('Individual merit functions:')
        for i in range(e.shape[2]):
            L+=np.sum(np.abs(e[cutx:-cuty,cutx:-cuty,i])**2)/e.shape[0]**2
        #    print(i,np.sum(np.abs(e[cutx:-cuty,cutx:-cuty,i])**2)/e.shape[0]**2)
    elif e.ndim==2:
        if cut!=None:
            L=np.sum(np.abs(e[cutx:-cuty,cutx:-cuty])**2)/e.shape[0]**2
        else:
            L=np.sum(np.abs(e[:,:])**2)
    return L


def derHj2(zj,p,norma):
    """
    Derivative of the OFT with respect to Zernike coefficient 'a_j' using
    the correlation definition of math_func2
    Ref: Paxman et al 1992, Eq. C14
    Inputs:
        zj: Zernike 'j' in Noll's basis
        p: generalized pupil function
    """
    K=p.shape[2]
    derivHj=np.zeros(p.shape,dtype='complex128')
    for i in range(K):
        C1=mf.corr(p[:,:,i],zj*p[:,:,i])
        C2=mf.corr(zj*p[:,:,i],p[:,:,i]) #C2=np.conj(np.flip(C1))
        corr1_corr2=C1-C2 #Substraction
        derivHj[:,:,i]=1j*(corr1_corr2)/norma[i]
    return derivHj

def derEj(derHk,Ok,Hk,Q,gamma):
    """
    Derivative with respect to a_j of the Merit function 'E'
    for phase diversity optimization
    Refs: 1) Lofdahl and Scharmer 1994, Eq.12
          2) Paxman et al 1992, Eq. C14
    Input:
        Q: 2D array returned by Qfactor
        Ok: 3D array with the FFTs of each of the PD images with different defoc.
        Hk: 3D array with the OTFs of the system for each defocus
        derHk: 3D array with the derivatives of the OTF for each defocus

    Output:
        2D complex array with the derivative of E for construction of the
        classical (Löfdahl & Scharmer) merit function
        or 3D complex array with the derivatives of each error function of the
        combination of pairs of the K images
    """
    #Derivate of Q
    derQ=-Q**3*np.sum(gamma*np.real(np.conj(Hk)*derHk),axis=2)

    #Array containing the derivative of each error function
    if Ok.ndim==3:
        K=Ok.shape[2]
        nk=int((K-1)*K/2)
        derivE=np.zeros((Ok.shape[0],Ok.shape[1],nk),dtype='complex128')
        for i in range(derivE.shape[2]):
            derivE[:,:,i]=derQ*arraymerit(Ok,Hk)[:,:,i]+Q*arraymerit(Ok,derHk)[:,:,i]
    elif Ok.ndim==2:
        derivE=derQ*arraymerit(Ok,Hk)+Q*arraymerit(Ok,derHk)
    return derivE


def derej(derE):
    """
    Inverse transform of the derivative of 'E' with respect to a_j.
    Input:
        derE: 2D complex array with the derivative of E for construction of the
            classical (Löfdahl & Scharmer) merit function
            or 3D complex array with the derivative of E for each possible
            combination of pairs of the K images
    Output:
        derej: IFFT of derE (2D array or 3D depending of the shape of E)
    """
    n=derE.shape[0]
    if derE.ndim==3: #Array with individual merit functions of K diverities
        dere=np.zeros(derE.shape)
        for i in range(derE.shape[2]):
            dere0=ifftshift(derE[:,:,i])
            dere0=n*n*ifft2(dere0).real
            dere[:,:,i]=dere0-np.mean(dere0)
    elif derE.ndim==2: #If dim=2, no multidimensional array, just a n x n image
        dere=ifftshift(derE)
        dere=n*n*ifft2(dere).real
        dere-=np.mean(dere)
    return dere

def Aij(derivEi,derivEj,cut=None):
    """
    Elements of the A matrix (Lofdahl and Scharmer 1994, Eq 15)
    """
    Aij=derivEi*np.conj(derivEj)
    cutx=cut
    cuty=cut

    if cut!=None:
        Aij=np.sum(Aij[cutx:-cuty,cutx:-cuty])/Aij[cutx:-cuty,cutx:-cuty].shape[0]**2
    else:
        Aij=np.sum(Aij[:])/Aij.shape[0]**2
    return Aij.real

def b_i(derivEi,E,cut=None):
    """
    Elements of the b vector (Lofdahl and Scharmer 1994, Eq 16)
    """
    bi=derivEi*np.conj(E)
    cutx=cut
    cuty=cut
    if cut!=None:
        bi=np.sum(bi[cutx:-cuty,cutx:-cuty])/bi[cutx:-cuty,cutx:-cuty].shape[0]**2
    else:
        bi=np.sum(bi[:])/bi.shape[0]**2

    return bi.real

def loop_opt(tol,Jmin,Jmax,w_cut,maxnorm,maxit,a0,a_d,RHO,THETA,ap,Ok,\
gamma,nuc,N,disp=True,cut=None,method='svd',ffolder='',K=2,inst='tumag'):
    """
    Function that carries out the iterative loop for optimizing the aberration
    coefficients in the Phase Diversity process. It employs the focused and
    defocused images of the object and returns an array with the optimum
    aberration coefficients.
    Inputs:
        tol: tolerance criterion to stop iterations (change in the norm of a)
        Jmin: lowest Noll's Zernike index to be fitted
        Jmax: upper limit of Noll's Zernike. The highest Zernike index fitted
            is Jmax-1.
        w_cut: SVD cutoff as a fraction of the maximum singular value
        maxnorm: max norm of a allowed in the iterative process
        maxit: maximum number of iterations allowed for the iterative process
        a0: initial guess for the Zernike indices (usually, array with zeroes)
        a_d: phase difference between images. It can be a scalar (classic
            two-images PD, or an array of length K if we use K images)
        RHO, THETA: output of 'sampling2' function
        ap: output of 'aperture' function
        Ok: output of 'prepare_PD'. Consists of an array with
            the FFT of the input PD images. The first index
            must correspond to the 'focused' image. The images should be
            arranged consistently with the phase differences of 'a_d'
        disp (True or False): to print some information on the iterative process
        cut: number of edge pixels on the images that are cropped
        method: optimization method ('svd' is highly recommended)
        gamma: output of "prepare_PD". Level of noise among images.
        nuc: critical frequency of the telescope (computed at the beginning
            of this module)
        N: size of each PD image (defined at the beginning of this module)
        ffolder: subfolder where the results will be saved (parent folder is 'txt')
        K: number of images employed (2 for focused-defocused classic PD)
    Output:
        a_iteration_Jmax_...:txt file with information on the merit function,
            norm of the change of Zernike coefficients and the Zernike coefficients
            at each iterative step. When subfielding is applied, only
            the information corresponding to the subfield that takes
            longer to be optimized is saved.
        a_optimized_Jmax_...: file with information of the optimized
            Zernike coefficients. If subfielding is applied, then
            one file per subfield is saved. Subfields are enumerated
            with one label (k)
    """
    print('gamma:',gamma)
    a=np.copy(a0)




    #Reference defocus for comparing purposes
    if inst=='tumag':
        wvl,fnum,_=tumag_params()
    d0=8*wvl*(fnum)**2 #Defocus shift to achieve 1 lambda peak-to-valley WF error
    a_d0=np.pi*d0/(8*np.sqrt(3)*wvl*fnum**2) #Defocus coefficient [rads]
    tic=time.time()
    it=0
    norm_delta=maxnorm

    #Initialization of A and b
    if Jmin==2: #To fit also individual tip/tilt terms
        A=np.zeros((2*(K-2)+Jmax-Jmin,2*(K-2)+Jmax-Jmin),dtype='float64')
        b=np.zeros((2*(K-2)+Jmax-Jmin,1),dtype='float64')
    else:
        A=np.zeros((Jmax-Jmin,Jmax-Jmin),dtype='float64')
        b=np.zeros((Jmax-Jmin,1),dtype='float64')


    #Initialization of the array to be saved in file
    L=0
    param=np.array([it,L,norm_delta])
    vfile=np.concatenate((param.reshape(len(param),1),a))
    file=np.copy(vfile)
    a_last=a[(Jmin-1):]+1
    a2=np.zeros((Jmax-1,1))

    #Planning for pyfftw
    if flag_pyfftw==1:
        pf=pyfftw.empty_aligned((RHO.shape[1],RHO.shape[1]),dtype='complex128')
        pd=pyfftw.empty_aligned((RHO.shape[1],RHO.shape[1]),dtype='complex128')
        zi_pf=pyfftw.empty_aligned((RHO.shape[1],RHO.shape[1]),dtype='complex128')
        zi_pd=pyfftw.empty_aligned((RHO.shape[1],RHO.shape[1]),dtype='complex128')
        zj_pf=pyfftw.empty_aligned((RHO.shape[1],RHO.shape[1]),dtype='complex128')
        zj_pd=pyfftw.empty_aligned((RHO.shape[1],RHO.shape[1]),dtype='complex128')

    if Jmin==2: #Activate tip/tilt correction in OTF
        tiptilt=True
    else:
        tiptilt=False
    while True:
        it+=1
        if disp is True:
            print('\n')
            print('Iteration',it)
        #Focused and defocused OTFs
        Hk,normhk=OTF(a,a_d,RHO,THETA,ap,tiptilt=tiptilt,norm=True,K=K) #Focused OTFs
        norma=normhk[0] #Norm for correlations

        #Q factor, merit function and focused and defocused pupils
        Q=Qfactor(Hk,gamma,nuc,N)
        noise_filt=filter_sch(Q,Ok,Hk,gamma,nuc,N)
        Ok_filt=np.zeros((N,N,K),dtype='complex128')
        pk=np.zeros((N,N,K),dtype='complex128')
        for i in range(K):
            Ok_filt[:,:,i]=noise_filt*Ok[:,:,i] #Filtered FFT
            if Jmin==2: #To fit also individual tip/tilt terms
                apup=select_tiptilt(a,i,K) #Zernike coefficients to compute pupil
                pk[:,:,i]=pupil(apup,a_d[i],RHO,THETA,ap) #Pupil
            else:
                pk[:,:,i]=pupil(a,a_d[i],RHO,THETA,ap) #Pupil

        E=meritE(Ok_filt,Hk,Q)
        e=merite(E)
        L_last=L
        L=meritl(e,cut=cut)
        derivL=np.abs(L-L_last)/L
        L2=meritL(E)

        #A matrix and b vector
        A_last=np.copy(A)
        b_last=np.copy(b)


        for i in range(Jmin,Jmax):
            zi=zk.zernikej_Noll(i,RHO,THETA) #Starting from i=1 (offset)
            derHk_i=derHj2(zi,pk,normhk)
            derivEi=derEj(derHk_i,Ok_filt,Hk,Q,gamma)
            derivei=derej(derivEi)
            if Jmin==2: #To fit also individual tip/tilt terms
                if i<4:
                    for ii in range(K-1):
                        kk=pair_index(ii,K)
                        b[2*ii+i-Jmin]=b_i(derivei[:,:,kk],e[:,:,kk],cut=cut)
                else:
                    b[2*(K-2)+i-Jmin]=b_i(derivei,e,cut=cut)
            else:
                b[i-Jmin]=b_i(derivei,e,cut=cut)

            for j in range(i,Jmax):
                if i==j:
                    derivEj=derivEi
                    derivej=derivei
                else:
                    zj=zk.zernikej_Noll(j,RHO,THETA)
                    derHk_j=derHj2(zj,pk,normhk)
                    derivEj=derEj(derHk_j,Ok_filt,Hk,Q,gamma)
                    derivej=derej(derivEj)
                if Jmin==2:
                    if j<4: #Choose only tip tilts corresponding to the same pair of images
                        for ii in range(K-1):
                            for jj in range(K-1):
                                ki=pair_index(ii,K)
                                kj=pair_index(jj,K)
                                if ki==kj:
                                    A[2*ii+i-Jmin,2*jj+j-Jmin]=Aij(derivei[:,:,ki],derivej[:,:,kj],cut=cut)
                                    A[2*jj+j-Jmin,2*ii+i-Jmin,]=A[2*ii+i-Jmin,2*jj+j-Jmin] #A is symmetric

                    elif i<4 and j>=4: #Compute matrix elements corresponding to same pair of images
                        for ii in range(K-1):
                            for jj in range(K-1):
                                ki=pair_index(ii,K)
                                kj=pair_index(jj,K)
                                if ki==kj:
                                    A[2*ii+i-Jmin,2*(K-2)+j-Jmin]=Aij(derivei[:,:,ki],derivej[:,:,kj],cut=cut)
                                    A[2*(K-2)+j-Jmin,2*ii+i-Jmin]=A[2*ii+i-Jmin,2*(K-2)+j-Jmin] #A is symmetric

                    else: #For other elements, sum the metrics of all pairs of images
                        A[2*(K-2)+i-Jmin,2*(K-2)+j-Jmin]=Aij(derivei,derivej,cut=cut)
                        A[2*(K-2)+j-Jmin,2*(K-2)+i-Jmin]=A[2*(K-2)+i-Jmin,2*(K-2)+j-Jmin] #A is symmetric
                else:
                    A[i-Jmin,j-Jmin]=Aij(derivei,derivej,cut=cut)
                    A[j-Jmin,i-Jmin]=A[i-Jmin,j-Jmin] #A is symmetric

        #plt.imshow(A,origin='lower',vmin=-1e-4,vmax=1e-4,cmap='seismic')
        #plt.colorbar()
        #plt.show()
        #plt.close()
        #quit()


        #SVD decomposition
        delta_a=mf.svd_solve(A,-b,w_cut,method=method,rms_limit=maxnorm)

        #Next iteration or break the loop
        norm_last=norm_delta
        norm_delta=np.linalg.norm(delta_a)
        if disp is True:
            print('Norm delta_a:',round(norm_delta,3))

        if norm_delta>maxnorm:
            print('\nOptimization stopped: Norm of delta_a >',maxnorm)
            print('Defocus:',round(a_d/a_d0,2),' lambda')
            print('Number of iterations:',it-1)
            break
        elif it>maxit:
            if disp is True:
                print('\nMaximun number of iterations reached:',maxit)
            break

        else:
            flag=0
            a_last=np.copy(a[(Jmin-1):])
            if Jmin==2: #To fit also individual tip/tilt terms
                a[(Jmin+1):]+=delta_a
            else:
                a[(Jmin-1):]+=delta_a

            param=np.array([it,L,norm_delta])
            vfile=np.concatenate((param.reshape(len(param),1),a))
            file=np.column_stack((file,vfile))
        if tol>norm_delta and it>2:
            break
        if delta_a.all()==0: #See svd_solve in math_func2 to check when this happens
            break

    if disp is True:
        if tol>norm_delta:
            if it<maxit:
                print('Iterations to reach convergence:',it)
            #print('Solution: a = \n',a)
        toc=time.time()
        if toc-tic<60:
            print('Elapsed time:', toc-tic, 'sec')
        else:
            print('Elapsed time:', (toc-tic)/60, 'min')
    #Save txt file
    flabel=['it','L','norm']
    ll=-1
    if Jmin==2: #To save individual tip/tilt terms
        for i in range(len(a)):
            if i>(2*K):
                flabel.append('a%g'%(i+3-2*K))
            elif i>0:
                if i%2==1:
                    ll+=1
                    flabel.append('tiltx%g'%ll)
                else:
                    flabel.append('tilty%g'%ll)
            elif i==0:
                flabel.append('offset')
    else:
        for i in range(len(a)):
            flabel.append('a%g'%(i+1))

    file=np.column_stack((flabel,file))
    try:
        os.mkdir('./txt/'+ffolder)
    except FileExistsError:
        print('./txt/'+ffolder+' already created')

    if isinstance(a_d, np.ndarray) is True:
        filename='./txt/'+ffolder+'a_iteration_Jmax_%g_K_%g.txt'%(Jmax,K)
    else:
        filename='./txt/'+ffolder+'a_iteration_Jmax_%g_K_%g.txt'%(Jmax,K)
    np.savetxt(filename,file,fmt='%s',delimiter='\t')
    return a



def minimization(Jmin,Jmax,a0,a_d,RHO,THETA,ap,Ok,\
gamma,nuc,N,disp=True,cut=None,ffolder='',K=2,jac=True):
    """
    Function that optimizes the aberration
    coefficients in the Phase Diversity process with a non linear minimization
    method
    WARNING: not prepared to include tip/tilt terms
    """


    def merit_func(a):
        Hk,normhk=OTF(a,a_d,RHO,THETA,ap,tiptilt=False,norm=True,K=K) #Focused OTFs
        norma=normhk[0] #Norm for correlations

        #Q factor, merit function and focused and defocused pupils
        Q=Qfactor(Hk,gamma,nuc,N)
        noise_filt=filter_sch(Q,Ok,Hk,gamma,nuc,N)
        Ok_filt=np.zeros((N,N,K),dtype='complex128')
        #pk=np.zeros((N,N,K),dtype='complex128')
        for i in range(K):
            Ok_filt[:,:,i]=noise_filt*Ok[:,:,i] #Filtered FFT
            #pk[:,:,i]=pupil(a,a_d[i],RHO,THETA,ap) #Pupil

        #Calculation of merit function
        E=meritE(Ok_filt,Hk,Q)
        e=merite(E)
        L=meritl(e,cut=cut)


        #Additive constant to the merit function (no change during optimization)
        #for i in range(K):
        #    Oshift=np.fft.fftshift(Ok[:,:,i]) #No filtering (Paxman, eq. 21)
        #    ok=N**2*np.real(np.fft.ifft2(Oshift))
        #    L-=meritl(ok,cut=cut)
        return L

    def jac_L(a):
        b=np.zeros((Jmax-1,1),dtype='float64')
        Hk,normhk=OTF(a,a_d,RHO,THETA,ap,tiptilt=False,norm=True,K=K) #Focused OTFs
        norma=normhk[0] #Norm for correlations

        #Q factor, merit function and focused and defocused pupils
        Q=Qfactor(Hk,gamma,nuc,N)
        noise_filt=filter_sch(Q,Ok,Hk,gamma,nuc,N)
        Ok_filt=np.zeros((N,N,K),dtype='complex128')
        pk=np.zeros((N,N,K),dtype='complex128')
        for i in range(K):
            Ok_filt[:,:,i]=noise_filt*Ok[:,:,i] #Filtered FFT
            pk[:,:,i]=pupil(a,a_d[i],RHO,THETA,ap) #Pupil
        E=meritE(Ok_filt,Hk,Q)
        e=merite(E)

        #Calculation of b
        for i in range(Jmin,Jmax):
            zi=zk.zernikej_Noll(i,RHO,THETA)
            derHk_i=derHj2(zi,pk,normhk)
            derivEi=derEj(derHk_i,Ok_filt,Hk,Q,gamma)
            derivei=derej(derivEi)
            b[i-1]=b_i(derivei,e,cut=cut)
        return 2*b


    meth='L-BFGS-B'
    #opt={'ftol':1e-9,'gtol':1e-8}#
    #opt={'ftol':1e-11,'gtol':1e-10}
    opt={'ftol':1e-15,'gtol':1e-14}
    if jac is True:
        minim=scipy.optimize.minimize(merit_func,a0,method=meth,jac=jac_L,options=opt)
    else:
        minim=scipy.optimize.minimize(merit_func,a0,method=meth,options=opt)
    print(minim)
    return np.array([minim.x]).T #To return an array consisting of 1 column

def minimization_jitter(Ok,gamma,nuc,N,cut=None):
    """
    Function that optimizes the jitter term for a set of two images,
    one free from jitter and another one affected by it.
    We do not take into account derivatives of the OTF.
    The image free from jitter must be the one at index 0.
    """
    def merit_func(sigma):
        #OTFs for the jitter term
        Hk_jitt,_=OTF_jitt(sigma,nuc,N,norm=None)
       
        #Q factor and noise filtering
        Q=Qfactor(Hk_jitt,gamma,nuc,N,reg1=0,reg2=1)
        noise_filt=filter_sch(Q,Ok,Hk_jitt,gamma,nuc,N)
        Ok_filt=np.zeros((N,N,2),dtype='complex128')
        for i in range(2):
            Ok_filt[:,:,i]=noise_filt*Ok[:,:,i] #Filtered FFT

        #Calculation of merit function
        E=meritE(Ok_filt,Hk_jitt,Q)
        e=merite(E)
        L=meritl(e,cut=cut)
        #print('Merit function:',L)
        #If L=0 because of an error computing Q
        #set it to one to avoid the program to think that the convergence succeeded
        if L==0:
            L=1 
        return L


    meth='Nelder-Mead'#'L-BFGS-B' #'Nelder-Mead'
    #opt={'ftol':1e-9,'gtol':1e-8}#
    #opt={'ftol':1e-11,'gtol':1e-10}
    opt={'ftol':1e-15,'gtol':1e-14}
    minim=scipy.optimize.minimize(merit_func,[0,0,0],method=meth,options=opt,
                                  bounds=((0,None),(0,None),(0,None)))
    print(minim)
    return np.array([minim.x]).T #To return an array consisting of 1 column

def read_image(file,ext,num_im=0,norma='yes'):
    """
    Function that opens an image and resizes it to fill the NxN detector
    """
    if norma=='yes':
        print('Read image:',file,'is normalized by its mean !!!!!')
    if ext=='.sav':
        from scipy.io.idl import readsav
        s = readsav(file+ext)
        I = s.imagen
        I = np.array(I,dtype='float')
        return I
    elif ext=='.fits':
        from astropy.io import fits
        I=fits.open(file+ext)
        data = I[0].data


        #If first dimension corresponds to focused and defocused images
        if data.ndim==3:
            if data.shape[1]==data.shape[2]:#Reorder axes and crop the image
                data=np.moveaxis(data,0,-1)
                N=data.shape[1]
                xmax=data.shape[0]
                xcen=int(xmax/2)
                l2=int(N/2)
                data=data[(xcen-l2):(xcen+l2),(xcen-l2):(xcen+l2),:]
        elif data.ndim==4:
            return data
        #If last dimension corresponds to focused and defocused images
        if data.shape[0]==data.shape[1]:
            if data.ndim==3: #If 'data' contains 3 dimensions
                if num_im==0: #If num_im=0, then get all images
                    data = data[:,:,:]
                else: #Else, get only the image at index num_im
                    data = data[:,:,num_im]
            elif data.ndim==2: #If 'data' containes 2 dimensions
                data=data[:,:]
            if norma=='yes': #Each image is normalized to its mean
                if data.ndim==3:
                    nima=data.shape[2]
                    for i in range(nima):
                        datai=data[:,:,i]
                        data[:,:,i]=datai/np.mean(datai)


                    data=data[0:data.shape[0],0:data.shape[1],:]
                    data = np.array(data,dtype='float64')

                    #try:
                    #    defoc=I[1].data
                    #    return data,defoc
                    #except:
                    return data
                elif data.ndim==2:
                    data=data/np.mean(data)
                    data=data[0:data.shape[0],0:data.shape[1]]
                    data = np.array(data,dtype='float64')
                    return data
            else:
                I = data
                I = np.array(I,dtype='float64')
                return I
        #If first dimension corresponds to focused and defocused images
        elif data.shape[0]==2:
            data1=data[0,:,:]
            I1 = data1/np.mean(data1[:])
            I1 = np.array(I1,dtype='float64')
            I1=I1[0:data.shape[1],0:data.shape[2]]
            data2=data[1,:,:]
            I2 = data2/np.mean(data2[:])
            I2 = np.array(I2,dtype='float64')
            I2=I2[0:data.shape[1],0:data.shape[2]]
            return I1, I2
        #If focused and defocused images are contained in a single image
        elif data.shape[1]==2*data.shape[0]:
            data1=data[:,:data.shape[0]]

            I1 = data1/np.mean(data1[:])
            I1 = np.array(I1,dtype='float64')
            data2=data[:,data.shape[0]:]
            I2 = data2/np.mean(data2[:])
            I2 = np.array(I2,dtype='float64')
            return I1, I2

    elif ext=='.png':
        I=im.open(file+ext)
        #I=I.convert('L')
        I=np.array(I,dtype='float')
        I=I[:,:,-1]
        if norma=='yes':
            I = I/np.mean(I)
        return I
    else:
        I=im.open(file+ext)
        I=I.convert('L') #Converts to black and white (not valid for np arrays)
        I=np.array(I,dtype='float') #Converts to numpy array for operation
                                    # purposes
        return I

def save_image(ima,exitfile,folder='txt'):
    """
    This function save a numpy array with shape (N,N)
    into a FITS file
    """
    from astropy.io import fits
    hdu = fits.PrimaryHDU(ima)
    hdul=fits.HDUList([hdu])
    hdul.writeto('./'+folder+'/'+exitfile+'.fits',overwrite=True)
    return

def gaussian_noise(SNR,image):
    row,col = image.shape
    mu  = np.mean(image)
    sigma= mu*SNR
    gauss = np.random.normal(0,sigma,(row,col))
    gauss = gauss.reshape(row,col)
    noisy = image + gauss
    return noisy

def apod(nx,ny,perc):
   """
   Apodization window of size nx x ny. The input parameter
   perc accounts for the percentage of the window that is apodized
   """
   nx = int(nx)
   ny = int(ny)

   wx = np.ones(nx)
   wy = np.ones(ny)

   nxw = int(perc*nx/100.)
   nyw = int(perc*ny/100.)

   wi_x = 0.5*(1.-np.cos(np.pi*np.arange(0,nxw,1)/nxw))
   wi_y = 0.5*(1.-np.cos(np.pi*np.arange(0,nxw,1)/nxw))
   wx[0:nxw] = wi_x
   wx[nx-nxw:nx]= wi_x[::-1]
   wy[0:nyw] = wi_y
   wy[ny-nyw:ny]= wi_x[::-1]

   win = np.outer(wx,wy)
   return win



def scanning(data,Lsiz=128,cut=29):
    """
    This function returns an array with the subpatches of the focused
    and defocused images
    """
    #Lsiz=int(N) #Subframe size
    lsiz=int(Lsiz-2*cut) #Size of the subframe that is not overlapped to others
    if np.floor(data.shape[0]/lsiz)>2:
        i_max=int(np.floor(data.shape[0]/lsiz)-1)
    else:#To perform a 2 x 2 subfielding
        i_max=int(np.floor(data.shape[0]/lsiz))
    i_vec=np.arange(0,i_max)

    #Overlapping of subframes
    kk=-1
    if data.ndim<3:
        Nima=1
        data_subpatch=np.zeros((int(i_max**2),Lsiz,Lsiz))
    else:    
        Nima=data.shape[2]
        data_subpatch=np.zeros((int(i_max**2),Lsiz,Lsiz,Nima))
    for i in i_vec:
        xc=int(Lsiz/2+lsiz*i)
        x0=int(xc-Lsiz/2)
        xf=int(xc+Lsiz/2)
        for j in i_vec:
            kk+=1
            yc=int(Lsiz/2+lsiz*j)
            y0=int(yc-Lsiz/2)
            yf=int(yc+Lsiz/2)
            if x0<0 or y0<0:
                print('Error in scanning.py: x0 or y0 cannot be negative')
                quit()
            elif xf>data.shape[0] or yf>data.shape[1]:
                print('Error in scanning.py: xf or yf larger than of size')
                print('xf:',xf)
                print('yf:',yf)
            if Nima==1:
                data_subpatch[kk,:,:]=data[x0:xf,y0:yf]    
            else:
                for n in range(Nima):
                    data_subpatch[kk,:,:,n]=data[x0:xf,y0:yf,n]
    return data_subpatch

def prepare_PD(ima,nuc,N,wind=True,kappa=100):
    """
    This function calculates gamma for each subpatch, apodizes the subpatch
    and calculates the Fourier transform of the focused and defocused images
    """
    #Initialization of arrays
    if ima.ndim==3:
        Nima=ima.shape[2]
    elif ima.ndim==2: #The array contains a single image
        ima=ima[..., np.newaxis] #To create a 3rd dummy dimension
        Nima=1    
    gamma=np.zeros(Nima)
    gamma[0]=1
    Ok=np.zeros((ima.shape[0],ima.shape[1],Nima),dtype='complex128')

    #Calculation of gamma before apodizing
    Of=mf.fourier2(ima[:,:,0])
    for i in range(1,Nima):
        #Fourier transforms of the images 
        Ok[:,:,i]=mf.fourier2(ima[:,:,i])
        gamma[i]=noise_power(Of,nuc,N)/noise_power(Ok[:,:,i],nuc,N)

    #Normalization to get mean 0 and apodization
    per_apod=10 #Percentage of apodization
    if wind==True:
        wind=apod(ima.shape[0],ima.shape[1],per_apod) #Apodization of subframes
    else:
        wind=np.ones((ima.shape[0],ima.shape[1]))

    #Apodization and FFT of focused image
    susf=np.sum(wind*ima[:,:,0])/np.sum(wind)
    of=(ima[:,:,0]-susf)*wind

    plt.close()
    #Of=mf.fourier2(of)
    Of=fft2(of)
    Of=Of/(N**2)

    #Apodization and FFTs of each of the K images
    for i in range(1,Nima):
        susi=np.sum(wind*ima[:,:,i])/np.sum(wind)
        imak=(ima[:,:,i]-susi)*wind

        #Fourier transforms of the images and of the PSFs
        #Ok[:,:,i]=mf.fourier2(imak)
        Ok[:,:,i]=fft2(imak)
        Ok[:,:,i]=Ok[:,:,i]/(N**2)

        #Compute and correct the shifts between images
        error,row_shift,col_shift,Gshift=sf.dftreg(Of,Ok[:,:,i],kappa)
        print('Residual shift:',row_shift,col_shift)

        Ok[:,:,i]=Gshift #We shift Ok[:,:,i]
        #Shift to center the FTTs
        Ok[:,:,i]=fftshift(Ok[:,:,i])

    Ok[:,:,0]=fftshift(Of)
    return Ok, gamma, wind, susf


def object_estimate(ima,a,a_d,wind=True,cobs=0,cut=29,low_f=0.2,tiptilt=False,
                    noise='default',reg1=0,reg2=1,inst='tumag'):
    """
    This function restores an image or an array of images employing a given
    set of Zernike coefficients.
    Inputs:
        ima: 2D array with the image to be restored or 3D array with PD
            images at different focus positions.
        a: array containing the set of Zernike coefficients in Noll's order    
        a_d: 0 if 'ima' is a 2D array or Numpy array with the list of
            focus positions if 'ima' is a 3D array
        wind: True or False. True to apodize the image an False otherwise.
        cobs: central obscuration of the telescope (0 if off-axis).
        cut: number of pixels to be excluded near the edges of the image (to
            compute the merit function).
        low_f: low limit of the noise filter, as defined in 'filter_sch'         
        tiptilt: True or False. If True, the first Zernike coefficients refer to
            the Tip/Tilt term of each of the images. False by Default because
            the Tip/Tilt terms usually vary a lot along the image and images are
            aligned with subpixel accuracy by a cross-correlation technique in
            'prepare_PD'
        noise: 'default' to be computed as filt_scharmer. Otherwise, this variable
            should contain a 2x2 array with the filter.
        inst: instrument for which we compute the parameters

    Output:
        object: restored image
        susf: offset defined in 'prepare_PD' 
        noise_filt: noise filter applied to deconvolve the images        
    """
    #Pupil sampling according to image size
    if inst=='tumag':
        wvl,fnum,Delta_x=tumag_params()
    N=ima.shape[0]
    nuc,R=compute_nuc(N,wvl,fnum,Delta_x)
    ap=aperture(N,R,cobs=cobs)
    RHO,THETA=sampling2(N,R)

    #Fourier transform images
    Ok, gamma, wind, susf=prepare_PD(ima,nuc,N,wind=wind)

    if isinstance(a_d, np.ndarray): #If a_d is an array...
        if a_d[0]==a_d[1]:#In this case, the 2nd image is a dummy image
            if a_d[0]==0:
                gamma=[1,0] #To account only for the 1st image

    #OTFs
    Hk,_=OTF(a,a_d,RHO,THETA,ap,norm=True,K=Ok.shape[2],tiptilt=tiptilt)

    
    #Restoration
    Q=Qfactor(Hk,gamma,nuc,N,reg1=reg1,reg2=reg2)
    
    if noise=='default':
        noise_filt=filter_sch(Q,Ok,Hk,gamma,nuc,N,low_f=low_f)
    else:
        noise_filt=noise

    Nima=Ok.shape[2]
    for i in range(0,Nima):
        #Filtering in Fourier domain
        Ok[:,:,i]=noise_filt*Ok[:,:,i]

    #Compute and print merit function
    if Nima>1:    
        E=meritE(Ok,Hk,Q)
        e=merite(E)
        L=meritl(e,cut=cut)#/(Nima*(Nima-1)/2)
        if a.any()!=0:
            print('Optimized merit function (K=%g):'%Nima,L)

    #Restoration
    O=Ffactor(Q,Ok,Hk,gamma)

    Oshift=np.fft.fftshift(O)
    o=np.fft.ifft2(Oshift)
    #o=np.fft.fftshift(o)
    o_restored=N**2*o.real
    object=o_restored+susf
    return object,susf,noise_filt



def object_estimate_jitter(ima,sigma,a,a_d,cobs=0,wind=True,low_f=0.2,
                    tiptilt=False,noise='default',reg1=0,reg2=1,inst='tumag'):
    """
    This function restores the OTF part of a jittered image.
    Inputs:
        ima: 2D array with the image to be restored or 3D array with PD
            images at different focus positions.
        sigma: rms of the jitter to be used when computing the OTF.
        a: numpy array with aberrations ordered following Noll's
        a_d: array or list with defocuses ([0,0] if we use only one image)
        wind: True or False. True to apodize the image an False otherwise.
        cobs: central obscuration of the telescope (0 if off-axis).
        low_f: low limit of the noise filter, as defined in 'filter_sch'         
        noise: 'default' to be computed as filt_scharmer. Otherwise, this variable
            should contain a 2x2 array with the filter.
        inst: instrument for which we compute the parameters
    Output:
        object: restored image
        susf: offset defined in 'prepare_PD' 
        noise_filt: noise filter applied to deconvolve the images        
    """
    #Sampling according to image size
    if inst=='tumag':
        wvl,fnum,Delta_x=tumag_params()
    N=ima.shape[0]
    nuc,R=compute_nuc(N,wvl,fnum,Delta_x)
    ap=aperture(N,R,cobs=cobs)
    RHO,THETA=sampling2(N,R)

    #Fourier transform images
    Ok, gamma, wind, susf=prepare_PD(ima,nuc,N,wind=wind)
    
    #Set gamma[1] to zero, as we do want to reconstruct only the first imag
    if isinstance(a_d, np.ndarray): #If a_d is an array...
        if a_d[0]==a_d[1]:#In this case, the 2nd image is a dummy image
            if a_d[0]==0:
                gamma=[1,0] #To account only for the 1st image
    
    #Jitter OTF
    Hk_jitt,_=OTF_jitt(sigma,nuc,N,norm=False)
    Hk_jitt[:,:,0]=Hk_jitt[:,:,1] #To employ the OTF of the jittered image

    #Aberration OTF
    Hk_aberr,_=OTF(a,a_d,RHO,THETA,ap,norm=True,K=Ok.shape[2],tiptilt=tiptilt)
 

    #Total OTF
    if np.all(a==0): #If all aberrations are zero
        Hk=Hk_jitt
    else:
        Hk=0*Hk_aberr
        for i in range(Ok.shape[2]):    
            Hk[:,:,i]=Hk_jitt[:,:,i]*Hk_aberr[:,:,i]
    
    #Restoration
    Q=Qfactor(Hk,gamma,nuc,N,reg1=reg1,reg2=reg2)
    
    if noise=='default':
        if np.all(sigma==0):
            noise_filt=filter_sch(Q,Ok,Hk,gamma,nuc,N,low_f=low_f)
        else:
            #If we have jitter, we reduce the cut-off frequency to avoid
            #spurious artifacts close to the Nyquist frequency
            nuc2=int(0.5*nuc)
            noise_filt=filter_sch(Q,Ok,Hk,gamma,nuc2,N,low_f=low_f)  
    else:
        noise_filt=noise

    Nima=Ok.shape[2]
    for i in range(0,Nima):
        #Filtering in Fourier domain
        Ok[:,:,i]=noise_filt*Ok[:,:,i]

    #Restoration
    O=Ffactor(Q,Ok,Hk,gamma)

    Oshift=np.fft.fftshift(O)
    o=np.fft.ifft2(Oshift)
    #o=np.fft.fftshift(o)
    o_restored=N**2*o.real
    object=o_restored+susf
    return object,susf,noise_filt

def import_aberr(file='IMaX_aberr.txt',Jmax_ab=33):
    aberr_file=file
    IMaX_aberr=np.loadtxt(aberr_file)
    a_aberr=IMaX_aberr[:(Jmax_ab-1),1] #Aberrations induced in the image [rads]
    norm_aberr=np.linalg.norm(a_aberr)
    return a_aberr,norm_aberr

def input_params(input):
    return Jmax,dmin,dmax,deltad,NSR,w_cut,a6,Jmax_ab,width

def zernike_from_txt(path):
        filename=path
        #Import txt files
        data=np.genfromtxt(filename,delimiter='\t',unpack=False,dtype=None,\
        encoding='utf-8')
        names=np.array(data[:,0],dtype='str')
        values=data[:,1]

        #Obtain values from imported data
        a1_ind=np.argwhere(names=='a1')[0][0]
        a=np.array(values[a1_ind:],dtype='float64')
        return a

def MTF_ideal(x):
    MTF=2/np.pi*(np.arccos(x)-x*np.sqrt(1-x**2))
    return MTF


def open_dark(path,header):
    I=fits.open(path)
    dark=I[0].data
    dark_header=I[0].header
    if dark_header['IMGDIRX'] == 'YES':
        printc("Flipping the dark: ", color=bcolors.FAIL)
        dark = np.fliplr(dark)

    PXBEG1 = int(header['PXBEG1']) - 1
    PXEND1 = int(header['PXEND1']) - 1
    PXBEG2 = int(header['PXBEG2']) - 1
    PXEND2 = int(header['PXEND2']) - 1
    # CHECK NACC
    acc = int(header['ACCACCUM']) * int(header['ACCCOLIT'])
    acc_dark = int(dark_header['ACCACCUM']) * int(dark_header['ACCCOLIT'])
    if acc != acc_dark:
        printc('WARNING - NACC NOT IDENTICAL DURING DARK CORRECTION', color=bcolors.FAIL)
        printc('DARK NACC ', acc_dark, ' DATA NACC ', acc, color=bcolors.FAIL)

    dark = dark[PXBEG2:PXEND2 + 1, PXBEG1:PXEND1 + 1]
    return dark

def open_flat(path,header):
    I=fits.open(path)
    flat=I[0].data


    PXBEG1 = int(header['PXBEG1']) - 1
    PXEND1 = int(header['PXEND1']) - 1
    PXBEG2 = int(header['PXBEG2']) - 1
    PXEND2 = int(header['PXEND2']) - 1

    if flat.ndim==3:
        flat=flat[0,:,:] #Flat at continuum

    flat = flat[PXBEG2:PXEND2 + 1, PXBEG1:PXEND1 + 1]
    return flat


def defocus_array(refocus,stp,Nima,foc_index,index_pos):
    """
    Returns array with defocuses in radians
    """
    if refocus=='PD_new':
        a_d=np.array([-2*stp,-stp,stp,2*stp,0])
    else:
        a_d=np.arange(0,Nima*stp,stp) #Array with defocuses
        a_d=a_d-a_d[foc_index] #Offset for focused image to have a 0 defocus
    a_d=a_d[index_pos] #To get only the defocuses we are interested
    a_d=a_d*np.pi/(np.sqrt(3)) #Defocus in radians
    return a_d

def initial_guess(K,Jmin,Jmax,guess='zeros'):
    """
    Returns initial guess for the optimization of the merit function
    in radian units

    guess:
        'zeros' -> all terms are zero
        'Fatima' -> Z4, Z9 and Z10 obtained by Fatima at STP-195 (0.33 au)
        'STP-195' -> terms obtained by me at STP-195 (0.33 au)
        'STP-136' -> terms obtained by me at STP-136 (0.51 au)
    """
    if Jmin==2:
        a0=np.zeros((2*K+Jmax-3,1)) #Offset + 2K tip/tilt terms + (Jmax-3) terms >=4
    else:
        a0=np.zeros((Jmax-1,1))
    if guess != 'zeros':
        if guess=='STP-136':
            a0[:12,0]=[0,0,0,0.0511,0.0199,0.0592,0.0046,-0.0047,-0.0247,\
            0.1025,0.0669,-0.0047]
        elif guess=='Fatima':
            a0[:12,0]=[0,0,0,0.3,0,0,0,0,0,-0.2,0.2,0]
        elif guess=='STP-195':
            Zguess=[0,0,0,0.2251,0.0034,0.1163,0.0158,0.0196,-0.0158,\
            0.2068,0.1546,0.0053,-0.0069,0.0131,0.0043,0.0043,-0.0039,\
            0.007,0.0053,0.0086,0.0005]
            a0[:(Jmax-1),0]=Zguess[:(Jmax-1)]
        elif guess=='STP-229/0.497 au':
            Zguess=[0,0,0,0.1757,-0.0048,0.0627,0.0076,-0.0222,-0.0424,\
            0.2653,0.1044,-0.0084,0.0013,-0.0247,0.0056,-0.0043,-0.0007,\
            -0.0064,0.0052,0.0175,0.0059]
            a0[:(Jmax-1),0]=Zguess[:(Jmax-1)]
        elif guess=='STP-228/0.393 au':
            Zguess=[0.0,0.0,0.0,0.1906,-0.0076,0.161,0.0292,0.0286,0.0057,\
            0.2577,0.1097,0.0412,0.0084,0.0286,0.0112,-0.011,0.0138,0.0069,\
            0.0101,-0.0462,0.011]
            a0[:(Jmax-1),0]=Zguess[:(Jmax-1)]
        elif guess=='STP-250':
            Zguess=[0,0,0,0.3833,0.0082,0.106,-0.002,-0.0094,-0.0302,0.2902,\
            0.1706,0.0091,-0.006,-0.0134,0.0129,0.0033,0.0077,-0.0073,\
            0.0033,0.0333,-0.0187]
            a0[:(Jmax-1),0]=Zguess[:(Jmax-1)]
        a0=2*np.pi*a0 #To convert into radians
    return a0

def radial_profile(data, center):
    """
    This function returns the radial profile of a
    2D image.
    """
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile

def power_radial(ima):
    """
    Radial power of a 2D image
    """
    power=np.abs(mf.fourier2(ima)/(ima.shape[0])**2)**2
    xp0=int(ima.shape[0]/2)
    power_radial=radial_profile(power,[xp0,xp0])
    return power_radial


def retrieve_aberr(N,R,k_max,Jmax,cobs,txtfolder):
    """
    Function that retrieves the average of the Zernike coeeficients saved in the
    txt files of a given directory.
    Inputs:
        k_max: Number of txt files employed to average the wavefront
        Jmax: Maximum index of the zernike polynomials
        cobs: Diameter of the central obscuration as a percentage of the aperture
        txtfolder: Path of the txt files
    """
    av=np.zeros((Jmax-1,k_max))
    for k in range(k_max):
        """
        Results after optimization files
        """
        filename=txtfolder+'/a_optimized_Jmax_%g_k_%g.txt'%(Jmax,k)
        #Import txt files
        data=np.genfromtxt(filename,delimiter='\t',unpack=False,dtype=None,\
        encoding='utf-8')
        names=np.array(data[:,0],dtype='str')
        values=data[:,1]

        #Obtain values from imported data
        a1_ind=np.argwhere(names=='a4')[0][0]
        a=np.array(values[a1_ind:],dtype='float64')
        a=np.concatenate((np.array([0,0,0]),a)) #First three are offset and tiptilt
        av[:,k]=a



    #Find optimizations for which a=0 and delete them
    index_ceros=np.argwhere(np.sum(np.abs(av), axis=0)==0)
    av=np.delete(av, index_ceros, axis=1)

    #Average zernike coefficients and wavefront
    a_aver=np.mean(av,axis=1)
    norm_aver=2*np.pi/np.linalg.norm(a_aver)
    print('WFE RMS:lambda/',np.round(norm_aver,2))
    return a_aver

def padding(ima):
    """
    Function that performs a symmetric padding of 10 % of the pixels of
    the image to reduce edge effects when reconstructing it.
    """
    size=ima.shape[0]
    pad_width = int(size*10/(100-10*2))
    if ima.ndim==3:
        ima_pad = np.zeros((size+pad_width*2,size+pad_width*2,ima.shape[-1]))
        for i in range(ima.shape[-1]):
            ima_pad[:,:,i]=np.pad(ima[:,:,i], pad_width=((pad_width, pad_width), (pad_width, pad_width)), mode='symmetric')  
    elif ima.ndim==2:
        ima_pad=np.pad(ima, pad_width=((pad_width, pad_width), (pad_width, pad_width)), mode='symmetric')
    return ima_pad

def restore_ima(ima,zernikes,pd=0,low_f=0.2,reg1=0.05,reg2=1,cobs=32.4):
    """
    This function restores a 2D or a 3D image using a given set of Zernike 
    coefficients and a modified Wiener filter that includes a 
    regularization in the denominator of the type reg1*(nu/nuc)**reg2
    Input:
        ima: 2D array of dimensions NxN or 3D array of dimensions
            NxNxNima, where Nima is the number of PD images
        zernikes: list or 1D array of Zernike coefficients ordered
            following Noll's 1979 rule
        pd: phase diversity between the images in case ima is a 3D array
        low_f: lower threshold of the optimum (Wiener's) filter
        reg1: 1st parameter of regularization term
        reg2: 2nd parameter of the regularization term 
        cobs: percentage of the central obscuration
    Output:
        ima_rest: restored image over the whole FOV of the input image
        noise_filt: noise filter employed for the restoration
    """
    #Image padding to reduce edge effects
    ima_pad=padding(ima)
    cut=int(0.1*ima_pad.shape[0]) #To later select the non-padded region

    #If we select only one image of the series
    ima_rest,_,noise_filt=object_estimate(ima_pad,zernikes,pd,
                                                wind=True,cobs=cobs,cut=cut,
                                                low_f=low_f,reg1=reg1,reg2=reg2)
    ima_rest=ima_rest[cut:-cut,cut:-cut]
    return ima_rest, noise_filt
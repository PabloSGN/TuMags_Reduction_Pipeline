"""
Same as math_func except for the definition of the correlation in 'corr',
which is changed to coincide with the one given by IDL
"""
import numpy as np
from matplotlib import pyplot as plt
import scipy
from scipy.fftpack import fftshift, ifftshift, fft2, ifft2
from tqdm import tqdm
import shift_func as sf

flag=0
#try:
#    import pyfftw
#    flag=1
#except:
#    print('pyfftw not installed on computer')
#    flag=0



def fourier2(f,s=None):
    """
    This function calculates the Direct Fast Fourier Transform of an
    array and shifts it to center its spectrum. Input must be a real numpy
    array
    Input:
        I: (real) 2D numpy array (image)
    Output:
        O: 2D numpy array with the Fourier transform of the input
    """
    #if flag==1:
    #    fft2=pyfftw.interfaces.numpy_fft.fft2
    #    pyfftw.interfaces.cache.enable() #Turn on cache for optimum performance

    #F=fftshift(fft2(f,s=s))
    F=fftshift(fft2(f))
    return F

def ifourier2(F,s=None):
    """
    This function calculates the Inverse Fast Fourier Transform of an
    array and shifts it to center its spectrum. Input must be a real numpy
    array
    Input:
        I: (real) 2D numpy array (image)
    Output:
        O: 2D numpy array with the Fourier transform of the input
    """
    #if flag==1:
    #    fft2=pyfftw.interfaces.numpy_fft.fft2
    #    pyfftw.interfaces.cache.enable() #Turn on cache for optimum performance
    #f=fftshift(ifft2(F,s=s))
    f=fftshift(ifft2(F))
    return f

def convfft(a,b,norm=None):
    """
    Parameters:
        f,g: Numpy vectors or 2D matrices
        norm:{None,True}, optional. 'True' for normalization purpose. Default
             is None
    Computes the convolution of two REAL vectors 'a' and 'b'. The order is
    important if normalization is set to be True!. It uses
    'rfft' Numpy function. The family of rfft functions is
    designed to operate on real inputs, and exploits the Hermitian symmetry
    of the Fourier Transform by computing
    only the positive frequency components, up to and including the Nyquist
    frequency. Thus, n input points produce n/2+1 complex output points.
    The inverses of this family assumes the same symmetry of its input,
    and for an output of n points uses n/2+1 input points.
    """
    if norm==True:
        c=rfft(a)*rfft(b/np.sum(np.abs(b)))
    else:
        c=rfft(a)*rfft(b)
    c=irfft(c)
    c=fftshift(c)
    return c

def autocorr(f):
    """
    This function returns the autocorrelation of a vector or a 2D matrix.
    Not normalized for a vector input.
    """
    if np.ndim(f)==1:
        F=fft(f)
        F=fftshift(F)
        power=(np.abs(F))**2
        c=ifft(power)
        c=fftshift(c)
    else:
        n=f.shape[1]
        #F=np.fft.fft2(f)
        #F=np.fft.fftshift(F)
        F=fourier2(f)
        power=n*n*(np.abs(F))**2
        c=ifourier2(power)
        #c=np.fft.ifft2(power)
        #c=np.fft.fftshift(c)
    return c

def corr(f,g,norma=False):
    """
    This function returns the correlation of two vector or 2D matrices f and g.
    It is important to notice that the order MATTERS in correlations, in contrast
    to in convolution. The normalization factor is chosen from Bonet
    "crosscorr_c.pro".
    Parameters:
        f,g: Numpy vectors or 2D matrices
    """

    n=f.shape[1]
    F=fft2(f)
    G=fft2(g)
    power=n*n*np.conj(F)*G #Normalized correlation
    c=ifft2(power)
    norma_corr=np.abs(c[0,0])
    #c=np.real(c) #This is only true for the auto-correlation
    c=ifftshift(c)


    if norma==True:
        return norma_corr,c
    else:
        return c

def svd_solve(A,b,w_cut,method='svd',rms_limit=2):
    """
    This function solves the system of equations Ax=b by calculating the
    inverse of A using the SVD method: x=A^(-1)*b; A^(-1)=V*S^(-1)*U'
    Inputs:
        A: 2D array of dimensions nxm (n>=m)
        b: 1D array of dimensions n
        w_cut: cut-off frequency for singular values (fraction of the maximum).
        Diagonal elements S^(-1) are zero for the positions of S where its
        value is less than w_cut.
        method: svd or lstsq (also based on svd and more efficient)
    """
    #Limiting parameters
    nsing_limit=6

    #Initizalization
    rms=rms_limit+1
    nsing=nsing_limit+1


    if method=='svd':
        #SVD method
        U,S,Vt=np.linalg.svd(A)
        sigma=w_cut*np.max(S)
        Sinv=np.where(S<sigma,0,(1/S))

        nsing=np.sum(Sinv>0)
        print('initial nsing:',nsing)
        k=nsing+1
        if (nsing+1)>nsing_limit:
            while rms>rms_limit and k>nsing_limit:
                k-=1 #k = nsing at first iteration
                Sinv2=0*Sinv
                Sinv2[:k]=Sinv[:k]
                Ainv=np.dot(np.transpose(Vt)*Sinv2,np.transpose(U))

                #Check that Ainv is correct at first iteration
                #Ainv2=np.linalg.pinv(A,rcond=w_cut,hermitian=False)

                delta_a=np.dot(Ainv,b)
                rms=np.linalg.norm(delta_a)

            print('final nsing:',k)
            print('effective w_cut:',round(S[k-1]/S[0],3))
        else:
            print('WARNING: nsing<nsing_limit in first iteration')
            print('NO OPTIMIZATION IS PERFORMED')
            delta_a=0*b
    elif method=='lstsq':
        w0=w_cut
        while rms>rms_limit and nsing>nsing_limit: #(Julián tiene puesto 'or')
            delta_a,resid,rank,S=np.linalg.lstsq(A,b,rcond=w0)
            rms=np.linalg.norm(delta_a)
            nsing=rank
            w0+=0.001
            #print('RMS;',rms,'nsing:',nsing)
        print('nsing:',rank)
        print('w_cut',round(w0-0.001,4))
    return delta_a

def pol2cart(rho, phi):
    """
    This function converts from polar coordinates to cartesian coordinates
    Arguments:
        rho: radial coordinate
        phi: polar angle
    """
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def cart2pol(x, y):
    """
    This function converts from cartesian coordinates to polar coordinates
    Arguments: x,y coordinates
    """
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def rgb2gray(rgb):
    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b
    return gray

def aprop(a):
    """
    Prints the basic properties of a numpy array
    """
    s=np.shape(a)
    d=np.ndim(a)
    si=np.size(a)
    its=a.itemsize
    ts=a.nbytes
    dt=a.dtype
    print('Shape ',s)
    print('Dim ',d)
    print('Size ',si)
    print('Item size ',its, ' bytes')
    print('Total size ',ts,' bytes')
    print('Data type ',dt)

def realign(im,Nseries,N,M,show_shift=False,return_shift=False):
    x=np.zeros(Nseries)
    y=np.zeros(Nseries)
    print("Realigning images...")
    for i in tqdm(range(1,Nseries)):
        AD=sf.AD(im[:,:,0],im[:,:,i],N,M)

        minim=np.unravel_index(np.argmin(AD, axis=None), AD.shape)

        #Parabolic simple interpolation
        shift=sf.parquad(AD,N)
        x0,y0=np.array(minim)-int(np.floor(N/2))
        x[i]=-(x0+shift[0])
        y[i]=-(y0+shift[1])


    if Nseries==2:
        im[:,:,i]=np.roll(im[:,:,i], int(round(x[i])), axis=0)
        im[:,:,i]=np.roll(im[:,:,i], int(round(y[i])), axis=1)
        print('Shift focused-defocused:',x[1],y[1])
        plt.imshow(AD)
        plt.show()
        plt.close()
    else:
        im[:,:,i]=np.roll(im[:,:,i], int(round(x[i])), axis=0)
        im[:,:,i]=np.roll(im[:,:,i], int(round(y[i])), axis=1)


    if show_shift == True:
        plt.plot(x,label='X shift')
        plt.plot(y,label='Y shift')
        plt.legend()
        plt.show()
        plt.close()

    if return_shift is False:
        return im
    if return_shift is True:
        return im,x[1],y[1]

import numpy as np
import time
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from scipy.fftpack import fftshift, ifftshift, fft2, ifft2

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

def dft_fjbm(F,G,kappa,dftshift,nr,nc,Nr,Nc,kernr,kernc):
    """
    Calculates the shift between a couple of images 'f' and 'g' with subpixel
    accuracy by calculating the IFT with the matrix multiplication tecnique.
    Shifts between images must be kept below 1.5 'dftshift' for the algorithm
    to work.
    Input:
        F,G: ffts of images 'f' and 'g' without applying any fftshift
        kappa: inverse of subpixel precision (kappa=20 > 0.005 pixel precision)
    Output:
    """
    #DFT by matrix multiplication
    M=F*np.conj(G)
    CC=kernr @ M @ kernc
    CCabs=np.abs(CC)
    ind = np.unravel_index(np.argmax(CCabs, axis=None), CCabs.shape)
    CCmax=CC[ind]
    rloc,cloc=ind-dftshift
    row_shift=-rloc/kappa
    col_shift=-cloc/kappa
    rg00=np.sum(np.abs(F)**2)
    rf00=np.sum(np.abs(G)**2)
    error=np.sqrt(1-np.abs(CCmax)**2/(rg00*rf00))
    Nc,Nr=np.meshgrid(Nc,Nr)

    Gshift=G*np.exp(1j*2*np.pi*(-row_shift*Nr/nr-col_shift*Nc/nc))
    return error,row_shift,col_shift,Gshift

def QD(f,g,N,M):
    """
    This function calculates the Quadratic Difference (Carmona et al 2014,
    System model of an Image Stabilization System) between the reference
    image 'f' and the shifted image 'g' in a NxN matrix.
    Input:
        f: reference image
        g: shifted image
        N: axis length of QD matrix
        M: axis length of 'f' and 'g' for calculation of the QD
    """
    if N%2==0:
        QD=np.zeros((N+1,N+1))
    else:
        QD=np.zeros((N,N))
    dim=f.shape
    x0=round(dim[0]/2)
    y0=round(dim[1]/2)
    N2=int(np.floor(N/2))
    kk=-1
    for k in range(-N2,N2+1):
        kk+=1
        ll=-1
        for l in range(-N2,N2+1):
            ll+=1
            QD[kk,ll]=np.sum((g[(x0+k-M):(x0+k+M+1),(y0+l-M):(y0+l+M+1)]\
            -f[(x0-M):(x0+M+1),(y0-M):(y0+M+1)])**2)
    return QD

def AD(f,g,N,M):
    """
    This function calculates the Absolute Difference (Carmona et al 2014,
    System model of an Image Stabilization System) between the reference
    image 'f' and the shifted image 'g' in a NxN matrix.
    Input:
        f: reference image
        g: shifted image
        N: axis length of QD matrix
        M: axis length of 'f' and 'g' for calculation of the QD
    """
    if N%2==0:
        AD=np.zeros((N+1,N+1))
    else:
        AD=np.zeros((N,N))
    dim=f.shape
    N2=int(np.floor(N/2))
    xi=round((dim[0]-M)/2)
    xf=round((dim[0]+M)/2)+1

    if xi<N2:
        xi=N2
    if xf>(dim[0]-N2):
        xf=dim[0]-N2
    kk=-1
    for k in range(-N2,N2+1):
        kk+=1
        ll=-1
        for l in range(-N2,N2+1):
            ll+=1
            AD[kk,ll]=np.sum(np.abs(g[(xi+k):(xf+k),(xi+l):(xf+l)]-\
            f[xi:xf,xi:xf]))
    return AD

def parquad(AD,N,method='simple'):
    """
    This functions calculates parabolic interpolation for the
    Absolute Difference based on the
    quadratic interpolation method (Carmona et al 2014)
    or on the minimum square method (Löfdahl, 2010) and returns the
    shift between images with subpixel accuracy
    """
    x0,y0=np.unravel_index(np.argmin(AD, axis=None), AD.shape)
    QD=AD**2 #QD is not the quadratic difference defined in QD function
    if method=='minimum_square':
        a20=2*QD[x0+1,y0+1]+2*QD[x0+1,y0]+2*QD[x0+1,y0-1]-4*QD[x0,y0+1]\
        -4*QD[x0,y0]-4*QD[x0,y0-1]+2*QD[x0-1,y0+1]+2*QD[x0-1,y0]+2*QD[x0-1,y0-1]
        a02=2*QD[x0+1,y0+1]-4*QD[x0+1,y0]+2*QD[x0+1,y0-1]+2*QD[x0,y0+1]\
        -4*QD[x0,y0]+2*QD[x0,y0-1]+2*QD[x0-1,y0+1]-4*QD[x0-1,y0]+2*QD[x0-1,y0-1]
        a11=3*QD[x0+1,y0+1]-3*QD[x0+1,y0-1]-3*QD[x0-1,y0+1]+3*QD[x0-1,y0-1]
        a10=2*QD[x0+1,y0+1]+2*QD[x0+1,y0]+2*QD[x0+1,y0-1]-2*QD[x0-1,y0+1]\
        -2*QD[x0-1,y0]-2*QD[x0-1,y0-1]
        a01=2*QD[x0+1,y0+1]-2*QD[x0+1,y0-1]+2*QD[x0,y0+1]-2*QD[x0,y0-1]\
        +2*QD[x0-1,y0+1]-2*QD[x0-1,y0-1]
    elif method=='simple':
        a10=(1/2)*(QD[x0+1,y0]-QD[x0-1,y0])
        a01=(1/2)*(QD[x0,y0+1]-QD[x0,y0-1])
        a20=(1/2)*(QD[x0+1,y0]-2*QD[x0,y0]+QD[x0-1,y0])
        a02=(1/2)*(QD[x0,y0+1]-2*QD[x0,y0]+QD[x0,y0-1])
        a11=(1/4)*(QD[x0+1,y0+1]-QD[x0-1,y0+1]-QD[x0+1,y0-1]+QD[x0-1,y0-1])
    shiftx=(a01*a11-2*a02*a10)/(4*a02*a20-a11**2)
    shifty=(a10*a11-2*a20*a01)/(4*a02*a20-a11**2)
    return shiftx,shifty

def kernel(x,method='cubic',extra_par=1):
    """
    This function calculates the Kernel employed for the interpolation of
    the image. It admits three cubic interpolation methods: cubic, catmull_rom
    and mitchell_netravali (approximation). The algorithms are based on
    "Principles of Digital Image Processing" (Burger 2009)
    Available methods: cubic, catmull_rom, mitchel_netravali,cubic_general,
    Lanczos
    """
    absx=np.abs(x)
    if method=='cubic':
        if absx<1:
            w=absx**3-2*absx**2+1
        elif absx<2:
            w=-absx**3+5*absx**2-8*absx+4
        else:
            w=0
    elif method=='catmull_rom':
        if absx<1:
            w=0.5*(3*absx**3-5*absx**2+2)
        elif absx<2:
            w=0.5*(-absx**3+5*absx**2-8*absx+4)
        else:
            w=0
    elif method=='mitchell_netravali':
        if absx<1:
            w=(1/18)*(21*absx**3-36*absx**2+16)
        elif absx<2:
            w=(1/18)*(-7*absx**3+36*absx**2-60*absx+32)
        else:
            w=0
    elif method=='cubic_general':
        a=extra_par
        if absx<1:
            w=(-a+2)*absx**3+(a-3)*absx**2+1
        elif absx<2:
            w=-a*absx**3+5*a*absx**2-8*a*absx+4*a
        else:
            w=0
    elif method=='Lanczos':
        n=extra_par
        if absx==0:
            w=np.sinc(x)
        elif absx<n:
            w=np.sinc(x)*np.sin(np.pi*x/n)/(np.pi*x/n)
        else:
            w=0
    return w

def interp2d(I,deltax,deltay,method='cubic',extra_par=1):
    """
    This function interpolates in two dimensions an image 'I' in order to
    translate it with subpixel accuracy. It employes the kernels defined
    in 'kernel' function.
    Available methods: cubic, catmull_rom, mitchell_netravali,cubic_general,
    Lanczos
    """

    #Shift in integer pixel units
    if deltax>=0:
        xpixel=np.floor(deltax)
        if deltay<0:
            ypixel=np.ceil(deltay)
        elif deltay>=0:
            ypixel=np.floor(deltay)
    elif deltax<0:
        xpixel=np.ceil(deltax)
        if deltay<0:
            ypixel=np.ceil(deltay)
        elif deltay>=0:
            ypixel=np.floor(deltay)
    deltax-=xpixel #Decimal part of X shift
    deltay-=ypixel #Decimal part of Y shift
    I=np.roll(I,(-int(xpixel),-int(ypixel)),axis=(0,1)) #Pixel shift


    #Subpixel shift
    Nx,Ny=I.shape
    Iest=np.zeros((Nx,Ny))
    if method=='Lanczos':
        n=extra_par
        imin=-n+1
        imax=n+1
    else:
        imin=-1
        imax=3
    for j in range(imin,imax):
        p=0
        for i in range(imin,imax):
            p+=np.roll(I,(-i,-j),axis=(0,1))*\
            kernel(deltax-i,method=method,extra_par=extra_par)
        Iest[:,:]+=kernel(deltay-j,method=method,extra_par=extra_par)*p
    return Iest

def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

def quadcell(f,g,R):
    sz=f.shape
    Nx=sz[0]
    Nx2=np.int(Nx/2)

    #Reference image
    A2=np.sum(f[:Nx2,Nx2:]) #Quadrant 1
    A1=np.sum(f[:Nx2,:Nx2]) #Quadrant 2
    A3=np.sum(f[Nx2:,:Nx2]) #Quadrant 3
    A4=np.sum(f[Nx2:,Nx2:]) #Quadrant 4
    deltay1=(A1+A2-A3-A4)/(4*R)
    deltax1=(A1+A3-A2-A4)/(4*R)

    #Shifted image
    A2=np.sum(g[:Nx2,Nx2:]) #Quadrant 1
    A1=np.sum(g[:Nx2,:Nx2]) #Quadrant 2
    A3=np.sum(g[Nx2:,:Nx2]) #Quadrant 3
    A4=np.sum(g[Nx2:,Nx2:]) #Quadrant 4
    deltay2=(A1+A2-A3-A4)/(4*R)
    deltax2=(A1+A3-A2-A4)/(4*R)

    #Relative shift
    deltax=deltax2-deltax1
    deltay=deltay2-deltay1
    return deltax,deltay

def hdmi_sensor(f,g,R,l,plot='no'):
    """
    Limb sensor with 4 photodiodes whose centers are separated 2*R each one
    Input:
        f: reference image
        g: displaced image
        R: radius of the sun (pixels)
        l: linear size of each of the photodiodes (pixels)
        plot (optional): plot sensors over reference image
    """
    sz=f.shape
    Nx=sz[0]
    Nx2=np.int(Nx/2)
    l2=np.int(l/2)
    c1=[Nx2-R,Nx2]#Center of photodiode 1
    c2=[Nx2,Nx2+R]#Center of photodiode 2
    c3=[Nx2+R,Nx2]#Center of photodiode 3
    c4=[Nx2,Nx2-R]#Center of photodiode 4

    #Reference image
    A1=np.sum(f[(c1[0]-l2):(c1[0]+l2+1),(c1[1]-l2):(c1[1]+l2+1)]) # +Y diode
    A2=np.sum(f[(c2[0]-l2):(c2[0]+l2+1),(c2[1]-l2):(c2[1]+l2+1)]) # +X diode
    A3=np.sum(f[(c3[0]-l2):(c3[0]+l2+1),(c3[1]-l2):(c3[1]+l2+1)]) # -Y diode
    A4=np.sum(f[(c4[0]-l2):(c4[0]+l2+1),(c4[1]-l2):(c4[1]+l2+1)]) # -X diode
    deltay1=(A1-A3)/(2*l)
    deltax1=(A2-A4)/(2*l)
    if plot=='yes':
        fig,ax=plt.subplots(1)
        ax.imshow(f)
        rect1 = patches.Rectangle(((c1[0]-l2),(c1[1]-l2)),l,l,\
        linewidth=1,edgecolor='r',facecolor='none')
        rect2 = patches.Rectangle(((c2[0]-l2),(c2[1]-l2)),l,l,\
        linewidth=1,edgecolor='g',facecolor='none')
        rect3 = patches.Rectangle(((c3[0]-l2),(c3[1]-l2)),l,l,\
        linewidth=1,edgecolor='b',facecolor='none')
        rect4 = patches.Rectangle(((c4[0]-l2),(c4[1]-l2)),l,l,\
        linewidth=1,edgecolor='w',facecolor='none')
        ax.add_patch(rect1)
        ax.add_patch(rect2)
        ax.add_patch(rect3)
        ax.add_patch(rect4)
        plt.show()
        quit()


    #Shifted image
    A1=np.sum(g[(c1[0]-l2):(c1[0]+l2+1),(c1[1]-l2):(c1[1]+l2+1)]) # +Y diode
    A2=np.sum(g[(c2[0]-l2):(c2[0]+l2+1),(c2[1]-l2):(c2[1]+l2+1)]) # +X diode
    A3=np.sum(g[(c3[0]-l2):(c3[0]+l2+1),(c3[1]-l2):(c3[1]+l2+1)]) # -Y diode
    A4=np.sum(g[(c4[0]-l2):(c4[0]+l2+1),(c4[1]-l2):(c4[1]+l2+1)]) # -X diode
    deltay2=(A1-A3)/(2*l)
    deltax2=(A2-A4)/(2*l)

    #Relative shift
    deltax=deltax2-deltax1
    deltay=deltay2-deltay1
    return deltax,deltay

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def realign_subpixel(ima,accu=0.05):
    """
    This function aligns a series of images with subpixel images using the Sicairos
    method.
    Input:
     ima: 3D array of the type (Nx,Ny,Nima). Last dimension corresponds to the
        index of the image through the series
     accu: accuracy of the alignment in pixel units
    Output: returns the aligned 3D array
    """
    kappa=1/accu #Kappa factor defined in Sicairos method (1/fraction of pixel)
    Gshift=fft2(ima[:,:,0]) #FFT of the first image of the series
    print('Re-aligning images ...')  
    ima_aligned=0*ima
    for j in range(ima.shape[-1]):
        F0=Gshift
        F_comp=fft2(ima[:,:,j])
        error,row_shift,col_shift,Gshift=dftreg(F0,F_comp,kappa)
        print('Shift of image',j,':',row_shift,col_shift)
        ima_aligned[:,:,j]=np.real(ifft2(Gshift))
    return ima_aligned
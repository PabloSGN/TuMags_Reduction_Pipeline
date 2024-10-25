import numpy as np
from scipy import special as sc
import os
dirname = os.path.dirname(__file__)
zernike_equiv = np.loadtxt(os.path.join(dirname, 'j_to_Noll.txt'), dtype='int')


def radialpol(m,n,rho):
    """
    This funcion calculates the radial polynomials of the Zernike polynomials
    Arguments:
        m,n: azimuthal and radial degrees of the Zernike polynomials
        rho: radial coordinate normalized to the unit (vector or array)
    """
    l=(n-np.abs(m))/2
    t=(n+np.abs(m))/2
    if np.size(rho)==1:
        r=0
    else:
        r=np.zeros(rho.shape)
    for s in range(int(l)+1):
        r+=(((-1)**s*sc.factorial(n-s))/(sc.factorial(s)*\
        sc.factorial(t-s)*sc.factorial(l-s)))*rho**(n-2*s)
    return r

def kroneckerDelta(m,n):
    """
    This function calculates the kroneker delta for indices m and n
    """
    if m==n:
        delta=1
    else:
        delta=0
    return delta

def zernike(m,n,rho,theta):
    """
    This function calculates the Zernike polinomial of degree m,n
    Arguments:
        m,n: azimuthal and radial degrees of the Zernike polynomials
        rho: radial coordinate normalized to the unit (vector or array)
        theta: polar angle (vector or array)
    """
    N=np.sqrt((2*(n+1))/(1+kroneckerDelta(m,0)))
    if m>0:
        Z=N*radialpol(m,n,rho)*np.cos(m*theta)
    elif m<0:
        Z=N*radialpol(m,n,rho)*np.sin(-m*theta);
    else:
        Z=N*radialpol(m,n,rho)

    Z=np.roll(Z, 1, axis=1)
    Z=np.roll(Z, 1, axis=0)
    return Z

def zernikej_Noll(j,rho,theta):
    """
    Returns the zernike polinomial using the equivalence between single
    indices (j) and athimuthal (m) and radial (n) degrees.
    Ref: Thibos, L. N. et al (2002). "Standards for reporting the optical
    aberrations of eyes"
    """
    m=zernike_equiv[j-1,1]
    n=zernike_equiv[j-1,2]
    zj=zernike(m,n,rho,theta)
    return zj

def write_Noll(j):
    """
    Writes in a txt file the zernike polinomial using the equivalence between single
    indices (j) and athimuthal (m) and radial (n) degrees.
    Ref: Thibos, L. N. et al (2002). "Standards for reporting the optical
         aberrations of eyes"
    """
    k=0
    n=-1
    jlist=np.arange(1,j+1)
    mlist=[]
    nlist=[]
    while k<j:
        n+=1
        mvec=np.arange(-n,n+1)
        mvec_bool=np.where((n+mvec)%2==0,True,False) #n+m must be even
        mvec=np.extract(mvec_bool,mvec)
        absmvec=np.abs(mvec)
        indices=np.argsort(absmvec)
        mvec=mvec[indices]
        mvec2=np.copy(mvec)
        for m in mvec:
            k+=1
            if k<=j:
                if k%2!=0:
                    if m<0:
                        mzj=m
                        mlist.append(mzj)
                        nlist.append(n)
                    else:
                        mzj=-m
                        mlist.append(mzj)
                        nlist.append(n)
                else:
                    if m>=0:
                        mzj=m
                        mlist.append(mzj)
                        nlist.append(n)
                    else:
                        mzj=-m
                        mlist.append(mzj)
                        nlist.append(n)
                mvec2=np.delete(mvec2,np.where(mvec2==m))
            else:break
    return jlist,mlist,nlist

def zernikej(j,rho,theta):
    """
    Mapping of the azimuthal and radial indices into a single index 'j' following
    OSA standark numeration
    Ref: Thibos, L. N. et al (2002). "Standards for reporting the optical
         aberrations of eyes"
    """
    k=-1
    n=-1
    while k<j:
        n=n+1
        for m in np.arange(-n-1,n+1):
            if (n+m)%2==0:
                k=k+1
                if k>=j:break
    zj=zernike(m,n,rho,theta)
    return zj

def zernikej_to_mn(Jmin,Jmax):
    """
    Shows in the command line the equivalence between single indices (j) and
    athimuthal (m) and radial (n) degrees.
    Ref: Thibos, L. N. et al (2002). "Standards for reporting the optical
         aberrations of eyes"
    """
    k=-1
    n=-1
    print('i','m','n','\n------')
    while k<Jmax:
        n=n+1
        for m in np.arange(-n-1,n+1):
            if (n+m)%2==0:
                k=k+1
                if k>=Jmin and k<=Jmax:
                    print(k,m,n)
                else:break

def zk_Noll_show(j):
    """
    Shows in the command line the equivalence between single indices (j) and
    athimuthal (m) and radial (n) degrees.
    Ref: Thibos, L. N. et al (2002). "Standards for reporting the optical
         aberrations of eyes"
    """
    k=0
    n=-1
    print('i','n','m','\n------')
    while k<j:
        n+=1
        mvec=np.arange(-n,n+1)
        mvec_bool=np.where((n+mvec)%2==0,True,False) #n+m must be even
        mvec=np.extract(mvec_bool,mvec)
        absmvec=np.abs(mvec)
        indices=np.argsort(absmvec)
        mvec=mvec[indices]
        mvec2=np.copy(mvec)
        for m in mvec:
            k+=1
            if k<=j:
                if k%2!=0:
                    if m<0:
                        print(k,n,m)
                    else:
                        print(k,n,-m)
                else:
                    if m>=0:
                        print(k,n,m)
                    else:
                        print(k,n,-m)
                mvec2=np.delete(mvec2,np.where(mvec2==m))
            else:break

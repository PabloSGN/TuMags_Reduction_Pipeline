"""
This module contains the functions to correct for crosstalk
through the method described in Jaeggli et al. 2022. 
(https://doi.org/10.3847/1538-4357/ac6506)

https://github.com/sajaeggli/adhoc_xtalk/blob/main/Jaeggli_etal_2022ApJ_AdHoc_Xtalk.ipynb

"""
import logging

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from .logutils import log_memory
from scipy.stats import linregress

def generate_squares(radius,divisions: int = 6):
    square_centers = np.linspace(-radius,radius,divisions+1,endpoint=True)[1:] - radius / divisions
    X,Y = np.meshgrid(square_centers, square_centers)
    return np.vstack([X.ravel(), Y.ravel()])

def evaluate_crosstalk(data, verbose=False, 
                       pthresh=0.05, png=False,
                       n_sigma=5, ctmethod="linfit",
                       region=[0,-1,0,-1]):
    
    def minimize_rms(x, y):
        norm = np.mean(y)
        def objective(residual): return np.std((x - residual * y) / norm)
        result = minimize(objective, 0.0, method='BFGS', 
                          options={'maxiter': 100, 'gtol': 1e-8, 'disp': verbose})
        return result.x[0]

    def pfit(x, y, n_sigma=5):
        # p = np.polyfit(x, y, deg=1)
        # residuals = y - np.polyval(p, x)
        # return np.polyfit(x[mask], y[mask], deg=1)
        slope, intercept, _,_,_ = linregress(x, y)
        residuals = y - (slope * x + intercept)
        mask = np.abs(residuals) < (n_sigma * np.std(residuals))
        return linregress(x[mask], y[mask])

    # Extract subregion and flatten
    area_of_interest = data[:, region[0]:region[1],region[2]:region[3]].reshape(data.shape[0], -1)
    y, q, u, v = area_of_interest[0], area_of_interest[1], area_of_interest[2], area_of_interest[3]

    #Apply threshold to V map corrected from offset. (as in jaeggli)
    V_mean_corr = area_of_interest[3]-np.mean(area_of_interest[3])
    pmap = np.abs(V_mean_corr)/area_of_interest[0]
    notpolar = pmap < pthresh

    # Compute slopes and intercepts
    if ctmethod == "linfit":
        slope_q, intercept_q, _, _, _ = pfit(y[notpolar], q[notpolar], n_sigma)
        slope_u, intercept_u, _, _, _ = pfit(y[notpolar], u[notpolar], n_sigma)
        slope_v, intercept_v, _, _, _ = pfit(y[notpolar], v[notpolar], n_sigma)
    elif ctmethod == 'rms':
        slope_q = minimize_rms(y[notpolar], q[notpolar])
        slope_u = minimize_rms(y[notpolar], u[notpolar])
        slope_v = minimize_rms(y[notpolar], v[notpolar])
        intercept_q = np.mean(q[notpolar]) - slope_q * np.mean(y[notpolar])
        intercept_u = np.mean(u[notpolar]) - slope_u * np.mean(y[notpolar])
        intercept_v = np.mean(v[notpolar]) - slope_v * np.mean(y[notpolar])
    else:
        raise ValueError("Invalid method. Use 'linfit' or 'rms'.")

    if verbose:
        fig, axes = plt.subplots(1, 3, figsize=(16, 6))
        titles = ['Q', 'U', 'V']
        datav = [(q, slope_q, intercept_q, notpolar),
                (u, slope_u, intercept_u, notpolar),
                (v, slope_v, intercept_v, notpolar)]

        for ax, (comp, slope, intercept, idx), title in zip(axes, datav, titles):
            fit_line = slope * y[idx] + intercept
            hb = ax.hexbin(y, comp, gridsize=50, cmap='inferno', mincnt=1)
            ax.plot(y[idx], fit_line, color='green', label='Fit')
            ax.set_title(f'Scatter Plot with Fit ({title})')
            ax.set_xlabel('I')
            ax.set_ylabel(title)
            ax.legend()
            fig.colorbar(hb, ax=ax, orientation='vertical', label='Density')

        if png:
            plt.savefig(f"{png}_fit.png")
            plt.close()
        else:
            plt.show()

    return intercept_q, intercept_u, intercept_v, slope_q, slope_u, slope_v

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

#Function that computes the merit function for a given Mueller matrix
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
    new_stokes = np.einsum('ij,abj->abi',iMM, np.squeeze(bs))
    Nwaves=new_stokes.shape[1] #Number of wavelenth samples
    
    # Minimization criteria
    out = np.abs(np.sum(new_stokes[:,:,0]*new_stokes[:,:,3],axis=1)) + \
          np.abs(np.sum(new_stokes[:,:,0]*new_stokes[:,:,2],axis=1)) + \
          np.abs(np.sum(new_stokes[:,:,0]*new_stokes[:,:,1],axis=1))

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
        if wvl=='all':      
            out = np.abs(np.sum(stI*stQ,axis=1)/kappa1)\
                +np.abs(np.sum(stI*stU,axis=1)/kappa2)\
                +np.abs(np.sum(stI*stV,axis=1)/kappa3)
        else:
            out = np.abs(stQ[:,wvl]*np.sum(stI,axis=1)/kappa1)\
                +np.abs(stU[:,wvl]*np.sum(stI,axis=1)/kappa2)\
                +np.abs(stV[:,wvl]*np.sum(stI,axis=1)/kappa3)    
    elif method=='all_wvls':
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
                       region=[0,-1,0,-1],
                       method='jaeggli',
                       last_wvl=None,plots=False,
                       roi = [0,-1,0,-1],
                       norm_wave = -1,
                       verbose = False,
                       ctmethod='linfit',
                       MM1a = None):
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
    logging.info(f"Starting cross-talk correction")
    log_memory("Before crosstalk")

    if norm is True:
        #Normalization of data
        norm_factor=np.median(data[roi[0]:roi[1],roi[2]:roi[3],norm_wave,0])
        data=data/norm_factor

    if method == 'standard':
        #loop in wavelengths
        data_corrected = np.copy(data)
        iq = np.zeros((data.shape[0]))
        iu = np.zeros((data.shape[0]))
        iv = np.zeros((data.shape[0]))
        sq = np.zeros((data.shape[0]))
        su = np.zeros((data.shape[0]))
        sv = np.zeros((data.shape[0]))
        if last_wvl != 0:
            input_data = np.reshape(np.einsum('lpij->pijl',data[:last_wvl]),(data.shape[1],data.shape[2],data.shape[3]*(data.shape[0]+last_wvl)))
            iq,iu,iv,sq,su,sv = evaluate_crosstalk(input_data,verbose=verbose,pthresh=pthresh,ctmethod=ctmethod,region=region)
            data_corrected[:,1,:,:] = data_corrected[:,1,:,:]  - sq*data_corrected[:,0,:,:]  - iq
            data_corrected[:,2,:,:] = data_corrected[:,2,:,:]  - su*data_corrected[:,0,:,:]  - iu
            data_corrected[:,3,:,:] = data_corrected[:,3,:,:]  - sv*data_corrected[:,0,:,:]  - iv
            
        else:
            for wvli in range(data.shape[0]):
                iq[wvli],iu[wvli],iv[wvli],sq[wvli],su[wvli],sv[wvli] = evaluate_crosstalk(data[wvli,:,:,:],verbose=verbose,pthresh=pthresh,ctmethod=ctmethod,region=region)
                data_corrected[wvli,1,:,:] = data_corrected[wvli,1,:,:]  - sq[wvli]*data_corrected[wvli,0,:,:]  - iq[wvli]
                data_corrected[wvli,2,:,:] = data_corrected[wvli,2,:,:]  - su[wvli]*data_corrected[wvli,0,:,:]  - iu[wvli]
                data_corrected[wvli,3,:,:] = data_corrected[wvli,3,:,:]  - sv[wvli]*data_corrected[wvli,0,:,:]  - iv[wvli]

        return data_corrected, (iq,iu,iv,sq,su,sv)
    
    if method == 'jaeggli' or method == 'all_wvls' or method == 'corr':

        #Reorder axis to convert into dimensions: [x,y,wavelength,stokes]
        data=np.moveaxis(data,0,-1)
        data=np.moveaxis(data,0,-1)
        Nwaves=data.shape[2]


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

        log_memory("Before weak data")

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
            plt.show(block=False)
            plt.close()        

        if method == 'jaeggli':

            if MM1a is not None:
                iMM1a = np.linalg.inv(MM1a)
                data_corrected =  np.einsum('ij,abcj->abci', iMM1a, full_data)
            else:
                #Minimize merit function
                result = minimize(fitfunc1, initial_guess, args=weak_region,
                    options={'maxiter': 1000, 'disp': verbose})
                # Apply correction for I<->QUV cross-talk
                MM1a = polmodel1(result.x[0],result.x[1], result.x[2])
                iMM1a = np.linalg.inv(MM1a)
                data_corrected =  np.einsum('ij,abcj->abci', iMM1a, full_data)

        elif method == 'all_wvls' or method == 'corr':

            MM1a=np.zeros((Nwaves,4,4))   
            data_corrected=full_data.copy()

            for wvli in range(Nwaves):
                logging.info(f"wave: {wvli+1}/{Nwaves}")

                fun=lambda x: fitfunc1(x, weak_region, wvl=wvli, method=method)
                    
                #Minimize merit function
                result = minimize(fun, initial_guess,
                    options={'maxiter': 1000, 'disp': verbose})


                # Apply correction for I<->QUV cross-talk
                MM1a[wvli,:,:] = polmodel1(result.x[0],result.x[1], result.x[2])
                iMM1a = np.linalg.inv(MM1a[wvli,:,:])
                data_inverted =  np.einsum('ij,abcj->abci', iMM1a, full_data)
                data_corrected[:,:,wvli,:] = data_inverted[:,:,wvli,:]

        #Move again axis to original positoin
        data_corrected=np.moveaxis(data_corrected,0,-1)
        data_corrected=np.moveaxis(data_corrected,0,-1)

        # log_memory("Ending cross-talk correction")

        return data_corrected, MM1a

def fit_mueller_matrix_2d(data,pthresh=0.02,norm=False,
                       divisions=14,region=[0,-1,0,-1],
                       method='standard',
                       last_wvl=None,plots=False,
                       verbose = False):

    logging.info(f"Starting 2D cross-talk correction")

    s = data.shape[-1]
    cx = s//2
    cy = s//2
    size = s//2
    size2 = size//divisions #half the size of the square

    ndiv = generate_squares(size,divisions = divisions)

    if verbose:
        fig, ax = plt.subplots(figsize=(8,8))
        val = data[0,0,cx,cy]
        im = ax.imshow(data[0,0,:,:],cmap='gray',clim=(val-val*1.5,val+val*1.5))
        for i in range(divisions**2):
            square = plt.Rectangle((ndiv[0,i] + cx - size2,ndiv[1,i] + cy - size2), size2 * 2, size2 * 2 , color='r', fill=False)
            ax.add_patch(square)
        plt.colorbar(im)
        plt.show()

    intercept = np.zeros((data.shape[0],3,divisions**2))
    slope = np.zeros((data.shape[0],3,divisions**2))
    mmatrix = np.zeros((data.shape[0],4,4,divisions**2))
    data_ct2D = np.copy(data)

    for i,loop in enumerate(range(divisions**2)):
        from_x, to_x = np.round(ndiv[0,i] + cx - size2).astype(int) , np.round(ndiv[0,i] + cx + size2).astype(int)
        from_y, to_y = np.round(ndiv[1,i] + cy - size2).astype(int) , np.round(ndiv[1,i] + cy + size2).astype(int)
        logging.info(f"loop {i} from {divisions**2}")

        data_ct2D[:,:,from_y:to_y,from_x:to_x], result = fit_mueller_matrix(
            data[:,:,from_y:to_y,from_x:to_x],
            pthresh=pthresh,
            norm=norm,
            region=region,
            method=method,
            last_wvl=last_wvl,
            plots=plots,
            verbose = verbose)
        if method == 'standard':
            intercept[:, 0, loop], intercept[:, 1, loop], intercept[:, 2, loop] = result[0], result[1], result[2]
            slope[:, 0, loop], slope[:, 1, loop], slope[:, 2, loop] = result[3], result[4], result[5]
        if method in ['jaeggli', 'all_wvls', 'corr']:
            mmatrix[:,:,:,loop] = result

    logging.info(f"finishing 2D cross-talk correction")

    if method == 'standard':
        return  data_ct2D, intercept, slope
    if method in ['jaeggli', 'all_wvls', 'corr']:
        return  data_ct2D, mmatrix, 0
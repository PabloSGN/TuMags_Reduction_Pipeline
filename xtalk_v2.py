"""
This module contains the functions to correct for crosstalk
through the method described in Jaeggli et al. 2022. 
(https://doi.org/10.3847/1538-4357/ac6506)

https://github.com/sajaeggli/adhoc_xtalk/blob/main/Jaeggli_etal_2022ApJ_AdHoc_Xtalk.ipynb

"""
import logging
import numpy as np 
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from scipy.optimize import minimize
from scipy.stats import linregress
from scipy.ndimage import sobel, laplace, gaussian_filter
from typing import List, Tuple, Dict, Any, Optional
from astropy.io.fits import Header
import csv

def _generate_squares(radius,divisions: int = 6):
    square_centers = np.linspace(-radius,radius,divisions+1,endpoint=True)[1:] - radius / divisions
    X,Y = np.meshgrid(square_centers, square_centers)
    return np.vstack([X.ravel(), Y.ravel()])


def evaluate_crosstalk(data, verbose=False, 
                       pthresh=0.05, png=False,
                       n_sigma=1, ctmethod="linfit",
                       region=[0,-1,0,-1],
                       pthresh_intensity=0,
                       xcomp='I',                   # añade 'dualI'
                       return_fit=False,
                       channels=('Q','U','V'),
                       # ===== NUEVO: intensidades dual beam =====
                       dualI=None,                  # None o array (2,H,W): (I_cam0, I_cam1)
                       dual_labels=('I0','I1'),     # etiquetas para cam0/cam1
                       # ===== MODO LOCAL CON DERIVADAS (sin Gauss) =====
                       use_local=False,          # igual que antes (no extendido a dual)
                       local_order=1,            # 1 -> (X, Ix, Iy); 2 -> añade Laplaciano
                       deriv_sigma=0.0,          # suavizado previo SOLO para derivadas (0 = off)
                       local_ridge_lambda=0.0,    # regularización ridge (0 = off)
                       fit_grad =  False
                       ):
    """
    Evalúa cross-talk en una sub-región usando:
      - Modelo lineal simple:  Y ≈ a + b * X                (xcomp='I' o 'V')
      - Modelo lineal dual:    Y ≈ a + b0*I0 + b1*I1        (xcomp='dualI', dualI=(I0,I1))
      - Modelo local (derivs): Y ≈ a + b*X + c*∂xX + d*∂yX [+ e*∇²X]   (como antes)
    """
    # ---------------- helpers ----------------
    def minimize_rms(source, target):
        """slope 'a' minimizando std((target - a*source)/norm); intercept por medias."""
        norm = np.mean(np.abs(target)) + 1e-12
        def objective(a):
            return np.std((target - a * source) / norm)
        result = minimize(objective, 0.0, method='BFGS',
                          options={'maxiter': 100, 'gtol': 1e-8, 'disp': verbose})
        return result.x[0]

    def linfit_clip(x, y, n_sigma=5):
        """Ajuste lineal (1 predictor) con 1 paso de sigma-clipping."""
        slope, intercept, *_ = linregress(x, y)
        resid = y - (slope*x + intercept)
        std = np.std(resid) + 1e-12
        mask = np.abs(resid) < (n_sigma * std)
        slope2, intercept2, *_ = linregress(x[mask], y[mask])
        return intercept2, slope2

    def linfit_clip_multi(A, y, n_sigma=5):
        """
        Ajuste lineal múltiple (>=1 predictores) con 1 paso de sigma-clipping.
        A tiene columna 0 todo 1's (intercepto).
        """
        # LS inicial
        theta, *_ = np.linalg.lstsq(A, y, rcond=None)
        resid = y - A @ theta
        std = np.std(resid) + 1e-12
        mask = np.abs(resid) < (n_sigma * std)
        if mask.sum() < max(3, A.shape[1] + 1):
            return theta  # no hay suficientes puntos tras clipping
        theta2, *_ = np.linalg.lstsq(A[mask], y[mask], rcond=None)
        return theta2

    def ridge_solve(A, y, lam=0.0):
        """Resuelve (A^T A + lam*L) theta = A^T y, sin penalizar intercepto."""
        AtA = A.T @ A
        Aty = A.T @ y
        if lam > 0:
            L = np.eye(AtA.shape[0])
            L[0,0] = 0.0  # no penalizar intercepto
            AtA = AtA + lam * L
        return np.linalg.solve(AtA, Aty)

    def flatten_mask(M2d, *arrs2d):
        outs = []
        idx = np.where(M2d)
        for A in arrs2d:
            outs.append(A[idx].ravel())
        return outs

    def build_derivs(X2d):
        """Derivadas (opcionalmente suavizado) para modelo local."""
        Xb = gaussian_filter(X2d, sigma=deriv_sigma) if (deriv_sigma and deriv_sigma > 0) else X2d
        Ix = sobel(Xb, axis=1)
        Iy = sobel(Xb, axis=0)
        if local_order >= 2:
            L = laplace(Xb)
        else:
            L = None
        return Xb, Ix, Iy, L

    def predictor_linear(a, b):
        def f(X2d):
            return a + b * X2d
        return f

    def predictor_local(a, b, c, d, e, order, deriv_sigma):
        def f(X2d):
            Xb = gaussian_filter(X2d, sigma=deriv_sigma) if (deriv_sigma and deriv_sigma > 0) else X2d
            Ix = sobel(Xb, axis=1)
            Iy = sobel(Xb, axis=0)
            yhat = a + b*X2d + c*Ix + d*Iy
            if order >= 2:
                yhat = yhat + e*laplace(Xb)
            return yhat
        return f

    def predictor_dual(a, b0, b1):
        def f(I0_2d, I1_2d):
            return a + b0*I0_2d + b1*I1_2d
        return f

    def predictor_solograd(c, d, e, order, deriv_sigma):
        def f(X2d):
            Xb = gaussian_filter(X2d, sigma=deriv_sigma) if (deriv_sigma and deriv_sigma > 0) else X2d
            Ix = sobel(Xb, axis=1)
            Iy = sobel(Xb, axis=0)
            # CORREGIDO: quitamos X2d ya que a=0, b=0 en el residuo
            yhat = c*Ix + d*Iy
            if order >= 2:
                yhat = yhat + e*laplace(Xb)
            return yhat
        return f

    # --------------- region & máscaras ---------------
    y1, y2, x1, x2 = region
    sub = data[:, y1:y2, x1:x2].astype(float)
    I2d, Q2d, U2d, V2d = sub[0], sub[1], sub[2], sub[3]

    if (xcomp == 'dualI') or (dualI is not None):
        if dualI is None:
            raise ValueError("xcomp='dualI' requiere dualI=(2,H,W) con (I_cam0, I_cam1).")
        I0 = dualI[0, y1:y2, x1:x2].astype(float)
        I1 = dualI[1, y1:y2, x1:x2].astype(float)

    eps = 1e-12
    Q_mc, U_mc, V_mc = Q2d - np.nanmean(Q2d), U2d - np.nanmean(U2d), V2d - np.nanmean(V2d)

    # máscara de selección (igual que antes)
    if np.isscalar(pthresh):
        if xcomp in ('I', 'dualI'):
            pmap = np.abs(V_mc) / (I2d + eps)
            sel2d = pmap < float(pthresh)
        elif xcomp == 'V':
            pmapV = np.abs(V_mc) / (I2d + eps)
            sel2d = pmapV >= float(pthresh)
        else:
            raise ValueError("xcomp debe ser 'I', 'V' o 'dualI'.")
    else:
        thr = np.asarray(pthresh, dtype=float)
        if thr.size != 3:
            raise ValueError("pthresh escalar o iterable de tres números (Q,U,V).")
        pmap_q = np.abs(Q_mc) / (I2d + eps)
        pmap_u = np.abs(U_mc) / (I2d + eps)
        pmap_v = np.abs(V_mc) / (I2d + eps)
        sel2d = np.ones_like(I2d, dtype=bool)
        if xcomp in ('I', 'dualI'):
            if thr[0] != 0: sel2d &= (pmap_q < thr[0])
            if thr[1] != 0: sel2d &= (pmap_u < thr[1])
            if thr[2] != 0: sel2d &= (pmap_v < thr[2])
        elif xcomp == 'V':
            if thr[0] != 0: sel2d &= (pmap_q >= thr[0])
            if thr[1] != 0: sel2d &= (pmap_u >= thr[1])
            if thr[2] != 0: sel2d &= (pmap_v >= thr[2])
        else:
            raise ValueError("xcomp debe ser 'I', 'V' o 'dualI'.")

    if pthresh_intensity != 0:
        sel2d &= (I2d > pthresh_intensity)

    # variable(s) predictoras
    if xcomp == 'I':
        X2d = I2d
        all_comps = {'Q': Q2d, 'U': U2d, 'V': V2d}
    elif xcomp == 'V':
        X2d = V2d
        all_comps = {'Q': Q2d, 'U': U2d}
    elif xcomp == 'dualI':
        all_comps = {'Q': Q2d, 'U': U2d, 'V': V2d}
    else:
        raise ValueError("xcomp debe ser 'I', 'V' o 'dualI'.")

    comps = {k: all_comps[k] for k in channels if k in all_comps}
    ch_list = list(comps.keys())

    # --------------- ajustes ---------------
    results = {}            # por canal: tuplas (a, b, c, d, e) o (a, b0, b1, 0, 0) en dual
    corrmap = {}            # por canal: {'I0': r, 'I1': r} si dualI
    for name in ch_list:
        Y2d = comps[name]

        # Modo local no lo extendemos a dual (mantiene X2d escalar)
        if not use_local and xcomp != 'dualI':
            # === Lineal (1 predictor) ===
            xv, yv = flatten_mask(sel2d, X2d, Y2d)
            if xv.size < 5:
                results[name] = (np.nan, np.nan, 0.0, 0.0, 0.0)  # a,b,c,d,e
                continue
            if ctmethod == 'linfit':
                a, b = linfit_clip(xv, yv, n_sigma)
            elif ctmethod == 'rms':
                b = minimize_rms(xv, yv); a = np.mean(yv) - b*np.mean(xv)
            else:
                raise ValueError("Método inválido. Usa 'linfit' o 'rms'.")
            c = d = e = 0.0
            results[name] = (float(a), float(b), float(c), float(d), float(e))

        elif not use_local and xcomp == 'dualI':
            # === Lineal múltiple (2 predictores: I0, I1) ===
            mask_ok = sel2d & np.isfinite(Y2d) & np.isfinite(I0) & np.isfinite(I1)
            if mask_ok.sum() < 10:
                results[name] = (np.nan, np.nan, np.nan, 0.0, 0.0)
                corrmap[name] = {dual_labels[0]: np.nan, dual_labels[1]: np.nan}
                continue
            yv = Y2d[mask_ok].ravel()
            i0 = I0[mask_ok].ravel()
            i1 = I1[mask_ok].ravel()

            # Matriz de diseño: [1, I0, I1]
            A = np.column_stack([np.ones_like(i0), i0, i1])

            if ctmethod == 'linfit':
                theta = linfit_clip_multi(A, yv, n_sigma=n_sigma)
            else:
                # Para dual, usamos LS incluso si ctmethod='rms'
                theta, *_ = np.linalg.lstsq(A, yv, rcond=None)

            a, b0, b1 = theta.tolist()
            results[name] = (float(a), float(b0), float(b1), 0.0, 0.0)

            # Correlaciones simples (Pearson) con cada haz
            def r_xy(x, y):
                x = x - np.nanmean(x); y = y - np.nanmean(y)
                den = (np.sqrt(np.nanmean(x*x)) * np.sqrt(np.nanmean(y*y)) + 1e-12)
                return float(np.nanmean(x*y) / den)
            corrmap[name] = {
                dual_labels[0]: r_xy(yv, i0),
                dual_labels[1]: r_xy(yv, i1)
            }
        elif use_local and not fit_grad:
            # === Local (derivadas, como antes; NO dual) ===
            if xcomp == 'dualI':
                raise NotImplementedError("Modelo local no implementado para dualI (usa lineal múltiple).")
            # Build feature maps
            Xb, Ix, Iy, L = build_derivs(X2d)
            cols = [np.ones_like(X2d), X2d, Ix, Iy]
            if local_order >= 2:
                cols.append(L)
            mats = flatten_mask(sel2d, *cols, Y2d)
            *Xs, yv = mats
            A = np.column_stack([x.ravel() for x in Xs])
            theta = ridge_solve(A, yv, lam=float(local_ridge_lambda))
            a, b, c, d = theta[:4]
            e = theta[4] if (local_order >= 2 and theta.size >= 5) else 0.0
            results[name] = (float(a), float(b), float(c), float(d), float(e))
        elif use_local and fit_grad:
            # === Local (derivadas sobre el residuo; NO dual) ===
            if xcomp == 'dualI':
                raise NotImplementedError("Modelo local no implementado para dualI (usa lineal múltiple).")
            # Build feature maps
            _, Ix, Iy, L = build_derivs(X2d)

            cols = [Ix, Iy]   # quitamos 1 (a) y X2d (b)
            if local_order >= 2:
                cols.append(L)

            mats = flatten_mask(sel2d, *cols, Y2d)
            *Xs, yv = mats
            A = np.column_stack([x.ravel() for x in Xs])
            theta = ridge_solve(A, yv, lam=float(local_ridge_lambda))

            # CORREGIDO: Forzamos a y b porque ajustamos sobre residuos
            a, b = 0.0, 1.0
            c, d = theta[:2]
            e = theta[2] if (local_order >= 2 and theta.size >= 3) else 0.0
            results[name] = (float(a), float(b), float(c), float(d), float(e))

        else:
            pass

    # --------------- tupla retrocompatible ---------------
    aQ,bQ,cQ,dQ,eQ = results.get('Q', (np.nan,)*5)
    aU,bU,cU,dU,eU = results.get('U', (np.nan,)*5)
    aV,bV,cV,dV,eV = results.get('V', (np.nan,)*5)

    if xcomp == 'dualI':
        # no hay única pendiente -> NaN en tuple para no inducir a error
        intercept_q, slope_q = aQ, np.nan
        intercept_u, slope_u = aU, np.nan
        intercept_v, slope_v = aV, np.nan
    else:
        intercept_q, slope_q = aQ, bQ
        intercept_u, slope_u = aU, bU
        intercept_v, slope_v = (aV, bV) if xcomp=='I' else (np.nan, np.nan)

    # --------------- diagnóstico (opcional) ---------------
    if verbose:
        import matplotlib.pyplot as plt
        n = len(ch_list)
        fig, axes = plt.subplots(1, n, figsize=(6*n, 5))
        if n == 1:
            axes = [axes]
        for ax, name in zip(axes, ch_list):
            a,b,c,d,e = results[name]
            if xcomp == 'dualI':
                # Y vs ŷ con dual
                I0p, I1p = I0, I1
                Y2d = comps[name]
                Yhat = a + b*I0p + c*I1p
                M = sel2d & np.isfinite(Y2d) & np.isfinite(Yhat)
                ax.hexbin(Y2d[M], Yhat[M], gridsize=60, cmap='inferno', mincnt=1)
                mn = np.nanmin(Y2d[M]); mx = np.nanmax(Y2d[M])
                ax.plot([mn, mx], [mn, mx], 'c--', lw=2, label='y = ŷ')
                rr0 = corrmap.get(name, {}).get(dual_labels[0], np.nan)
                rr1 = corrmap.get(name, {}).get(dual_labels[1], np.nan)
                ax.set_xlabel(f'{name} [DN]'); ax.set_ylabel(f'ŷ_{name} [DN]')
                ax.set_title(f'{name}: dualI (r_{dual_labels[0]}={rr0:.3f}, r_{dual_labels[1]}={rr1:.3f})')
                ax.legend()
            elif not use_local:
                # lineal simple (como antes)
                xv, yv = flatten_mask(sel2d, X2d, comps[name])
                xx = np.linspace(np.nanmin(xv), np.nanmax(xv), 200)
                ax.hexbin(xv, yv, gridsize=60, cmap='inferno', mincnt=1)
                ax.plot(xx, a + b*xx, 'c-', lw=2, label=f'{name} = a + b·{xcomp}')
                ax.set_xlabel(f'{xcomp} [DN]'); ax.set_ylabel(f'{name} [DN]')
                ax.set_title(f'{name} vs {xcomp} (linear)')
                ax.legend()
            else:
                # local: Y vs ŷ
                Xb, Ix, Iy, L = build_derivs(X2d)
                Y2d = comps[name]
                Yhat = a + b*X2d + d*Iy + c*Ix + (e*L if local_order >= 2 else 0.0)
                M = sel2d & np.isfinite(Y2d) & np.isfinite(Yhat)
                ax.hexbin(Y2d[M], Yhat[M], gridsize=60, cmap='inferno', mincnt=1)
                mn = np.nanmin(Y2d[M]); mx = np.nanmax(Y2d[M])
                ax.plot([mn, mx], [mn, mx], 'c--', lw=2, label='y = ŷ')
                ax.set_xlabel(f'{name} [DN]'); ax.set_ylabel(f'ŷ_{name} [DN]')
                ax.set_title(f'{name}: local (order={local_order})')
                ax.legend()

        plt.tight_layout()
        if isinstance(png, str) and len(png) > 0:
            plt.savefig(f"{png}_fit_{xcomp}.png", dpi=150); plt.close()
        else:
            plt.show()

    if return_fit:
        fit_dict = {}
        for name, (a,b,c,d,e) in results.items():
            if xcomp == 'dualI':
                pred_ct = predictor_dual(a, b, c)         # b->b0, c->b1 en este storage
                def correct_dual(Y2d, I0_2d, I1_2d):
                    return Y2d - pred_ct(I0_2d, I1_2d)
                fit_dict[name] = {
                    'a': a, 'b0': b, 'b1': c,             # ojo: aquí c es b1
                    'predict_ct_dual': pred_ct,
                    'correct_dual': correct_dual,
                    'xcomp': 'dualI',
                    'model': 'linear',
                    'deriv_sigma': float(deriv_sigma),
                    'lambda': float(local_ridge_lambda),
                    'corr': corrmap.get(name, {})
                }
            elif not use_local and not fit_grad:
                pred = predictor_linear(a, b)
                fit_dict[name] = {
                    'a': a, 'b': b,
                    'predict': pred,
                    'xcomp': xcomp,
                    'model': 'linear',
                    'deriv_sigma': float(deriv_sigma),
                    'lambda': float(local_ridge_lambda)
                }
            elif use_local and not fit_grad:
                pred = predictor_local(a, b, c, d, e, order=local_order, deriv_sigma=deriv_sigma)
                fit_dict[name] = {
                    'a': a, 'b': b, 'c': c, 'd': d, 'e': e,
                    'predict': pred,
                    'xcomp': xcomp,
                    'model': f'local{local_order}',
                    'deriv_sigma': float(deriv_sigma),
                    'lambda': float(local_ridge_lambda)
                }
            elif use_local and fit_grad:
                pred = predictor_solograd(c, d, e, order=local_order, deriv_sigma=deriv_sigma)
                fit_dict[name] = {
                    'a': a, 'b': b, 'c': c, 'd': d, 'e': e,
                    'predict': pred,
                    'xcomp': xcomp,
                    'model': f'grad{local_order}',
                    'deriv_sigma': float(deriv_sigma),
                    'lambda': float(local_ridge_lambda)
                }
            else:
                pass
        return (intercept_q, intercept_u, intercept_v, slope_q, slope_u, slope_v), fit_dict

    return intercept_q, intercept_u, intercept_v, slope_q, slope_u, slope_v
        
# JAEGGLI STUFF
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
def  minimize_for_model1(iMM,bs,wvl='all',method='jaeggli'):
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

def fit_mueller_matrix(
    data: np.ndarray,
    strategy: str = "simultaneous",
    pthresh: float = 0.02,
    norm: bool = False,
    region: List[int] = [0, -1, 0, -1],     # (y1,y2,x1,x2)
    method: str = 'standard',               # 'standard' | 'jaeggli' | 'all_wvls' | 'corr'
    last_wvl: Optional[int] = None,         # None -> todas las λ; int -> usa data[:last_wvl]
    plots: bool = False,
    roi: List[int] = [0, -1, 0, -1],        # (y1,y2,x1,x2) para normalización
    norm_wave: int = -1,                    # índice espectral para normalizar
    verbose: bool = False,
    ctmethod: str = 'linfit',               # 'linfit' | 'rms' (se pasa a evaluate_crosstalk)
    MM1a: Optional[np.ndarray] = None,
    # ===== NUEVO: intensidades dual beam =====
    dualI=None,                  # None o array (2,H,W): (I_cam0, I_cam1)
    # --- opciones específicas del flujo STANDARD (crosstalk con evaluate_crosstalk) ---
    use_local: bool = False,                # True -> modelo con derivadas
    # fit_grad: bool = False,                 # True -> fuerza a=0, b=0 para ajustar solo el gradiente (residuo)
    local_order: int = 1,                   # 1 -> (Ix,Iy), 2 -> (+ Laplaciano)
    deriv_sigma: float = 0.0,               # suavizado previo de derivadas
    local_ridge_lambda: float = 0.0,        # regularización para modelo local
    aggregate_wavelengths: bool = False,    # True -> ajusta una sola vez apilando todas las λ en 2D
    channels: Tuple[str, ...] = ('Q','U','V'),  # qué canales corregir en STANDARD
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Fit and correct cross-talk in spectro-polarimetric data.

    Parameters
    ----------
    data : np.ndarray
        Array with Stokes parameters AFTER dual-beam demodulation.
        Shape: (wavelength, stokes, x, y) with stokes ordered as [I,Q,U,V].
    pthresh : float
        Threshold on fractional polarization (Jaeggli-like criterion) used by evaluate_crosstalk.
    norm : bool
        If True, normalize data by median(I) in ROI at 'norm_wave'.
    region : [y1, y2, x1, x2]
        Image subregion for cross-talk estimation (and/or merit minimization).
        Negative indices follow Python slicing (e.g., -1 means end).
    method : {'standard', 'jaeggli', 'all_wvls', 'corr'}
        - 'standard': call evaluate_crosstalk (linear or local) and apply correction.
        - 'jaeggli' / 'all_wvls' / 'corr': fit diattenuation Mueller matrix (requires fitfunc1, polmodel1).
    last_wvl : int or None
        If not None, only wavelengths [:last_wvl] are used for estimation (Jaeggli path)
        or for aggregation (if aggregate_wavelengths=True).
    plots : bool
        If True, produce diagnostic plots in Jaeggli path (unchanged from your code).
    roi : [y1, y2, x1, x2]
        ROI for normalization when norm=True.
    norm_wave : int
        Spectral index used to compute normalization factor (median I over ROI).
    verbose : bool
        Verbosity flag (passed to evaluate_crosstalk / optimizer).
    ctmethod : {'linfit','rms'}
        Fitting method passed to evaluate_crosstalk (standard path).
    MM1a : np.ndarray or None
        If provided in 'jaeggli' method, it will be used (no optimization).
    use_local : bool
        STANDARD path: use local model with derivatives (no Gaussian) in evaluate_crosstalk.
    local_order : {1,2}
        STANDARD path: 1 -> (Ix,Iy), 2 -> adds Laplacian term.
    deriv_sigma : float
        STANDARD path: Gaussian smoothing of X for derivatives only.
    local_ridge_lambda : float
        STANDARD path: Ridge regularization for local model (no penalty on intercept).
    aggregate_wavelengths : bool
        STANDARD path: if True, stack all wavelengths (or [:last_wvl]) into a single 2D plane
        and fit once (global coefficients). If False, fit per wavelength independently.
    channels : tuple of {'Q','U','V'}
        STANDARD path: which channels to estimate and correct.

    Returns
    -------
    data_corrected : np.ndarray
        Corrected Stokes cube with same shape as 'data' (wavelength, stokes, x, y).
    info : dict
        Metadata and results:
            - method, region, normalized, normalization_factor
            - (STANDARD)
                * mode: 'global' or 'per_wavelength'
                * coeffs_global: dict or None
                * coeffs_per_wvl: list[dict] (one per λ) or None
                * pixel_counts: int or list[int]
            - (JAEGGLI/ALL_WVLS/CORR)
                * MM1a: np.ndarray (matrix or per-wavelength)
                * optimizer_result: optional (not stored here by default)
    """

    # -------------------- Validaciones básicas --------------------
    if data.ndim != 4 or data.shape[1] != 4:
        raise ValueError("Expected data shape (wavelength, 4, x, y) with Stokes=[I,Q,U,V].")

    def _normalize_last_wvl(last_wvl, n_wvl):
        """
        Normaliza last_wvl a:
        - None  -> usa todas las longitudes
        - int k -> usa [:k] (k > 0)
        Acepta entradas tipo: None, '', 'None', 'none', 'all', -1, 0, enteros.
        Si k excede n_wvl, lo recorta a n_wvl. Si k <=0, devuelve None.
        """
        if last_wvl is None:
            return None
        if isinstance(last_wvl, str):
            s = last_wvl.strip().lower()
            if s in ('', 'none', 'all', 'na', 'null'):
                return None
            try:
                v = int(s)
            except ValueError:
                return None
        else:
            v = int(last_wvl)

        if v <= 0:
            return None
        if v > n_wvl:
            return n_wvl
        return v

    def _norm_if_needed(arr: np.ndarray) -> Tuple[np.ndarray, float]:
        if not norm:
            return arr, 1.0
        y1, y2, x1, x2 = roi
        I0 = arr[norm_wave, 0, y1:y2, x1:x2]
        nf = float(np.nanmedian(I0))
        if nf == 0 or not np.isfinite(nf):
            logging.warning("Normalization factor is zero or non-finite; skipping normalization.")
            return arr, 1.0
        return arr / nf, nf

    def _slice_region(y1y2x1x2: List[int]) -> Tuple[slice, slice]:
        y1, y2, x1, x2 = y1y2x1x2
        return slice(y1, y2), slice(x1, x2)

    def _stack_wavelengths_for_ct(arr: np.ndarray, wvl_slice: slice, reg: List[int]) -> np.ndarray:
        """
        Convierte (λ,4,x,y) -> (4, X, Y*L) apilando λ a lo ancho.
        Es ideal para ajustar una sola vez con evaluate_crosstalk.
        """
        Ysl, Xsl = _slice_region(reg)
        sub = arr[wvl_slice, :, Ysl, Xsl]      # (L, 4, h, w)
        L, _, h, w = sub.shape
        plane = sub.transpose(1, 2, 3, 0).reshape(4, h, w * L)  # (4, h, w*L)
        return plane

    n_wvl = data.shape[0]
    last_wvl = _normalize_last_wvl(last_wvl, n_wvl)

    # -------------------- Normalización opcional --------------------
    data = data.astype(float, copy=False)
    data, norm_factor = _norm_if_needed(data)

    # -------------------- Región y recorte espectral --------------------
    reg = region[:]  # [y1,y2,x1,x2]
    ysl, xsl = _slice_region(reg)
    wsl = slice(None) if last_wvl is None else slice(0, last_wvl)

    # -------------------- Salidas comunes --------------------
    data_corrected = np.copy(data)
    info: Dict[str, Any] = {
        'method': method,
        'region': tuple(reg),
        'normalized': bool(norm),
        'normalization_factor': float(norm_factor),
    }

    if strategy == 'simultaneous':
        # fit conjunto Q,U,V
        # ======================== MÉTODO STANDARD =========================
        if method.lower() == 'standard':
            ec_kwargs = dict(
                xcomp='I',
                ctmethod=ctmethod,
                n_sigma=3,
                pthresh=pthresh,
                pthresh_intensity=0,
                region=[0, -1, 0, -1],  
                channels=tuple(channels),
                return_fit=True,
                fit_grad=False,
                use_local=use_local,
                local_order=local_order,
                deriv_sigma=deriv_sigma,
                local_ridge_lambda=local_ridge_lambda,
                verbose=verbose,
            )

            coeffs_per_wvl: Optional[List[Dict[str, Any]]] = None
            coeffs_global: Optional[Dict[str, Any]] = None
            pixel_counts: Any = None

            # ---------- Ajuste GLOBAL (apilando λ) ----------
            if aggregate_wavelengths:
                logging.warning("aggregate_wavelengths is True.")

                plane = _stack_wavelengths_for_ct(data[:, :, ysl, xsl], wsl, [0, -1, 0, -1])  # (4, h, w*L)
                (_, fit_glob) = evaluate_crosstalk(plane, **ec_kwargs)
                coeffs_global = fit_glob

                max_w = n_wvl if last_wvl is None else int(last_wvl)
                for i in range(max_w):
                    I_full = data[i, 0, :, :] 
                    for ch, idx in (('Q', 1), ('U', 2), ('V', 3)):
                        if ch not in channels or ch not in fit_glob:
                            continue
                        Y_full = data[i, idx, :, :]
                        Yhat_full = fit_glob[ch]['predict'](I_full)
                        M = np.isfinite(Y_full) & np.isfinite(Yhat_full)
                        out = data_corrected[i, idx, :, :]
                        out[M] = Y_full[M] - Yhat_full[M]

                info.update({
                    'mode': 'global',
                    'coeffs_global': coeffs_global,
                    'coeffs_per_wvl': None,
                    'pixel_counts': int(np.isfinite(plane[0]).sum())
                })
                return data_corrected, info

            # ---------- Ajuste POR LONGITUD DE ONDA ----------
            coeffs_per_wvl = []
            pixel_counts = []
            for i in range(n_wvl):
                if last_wvl is not None and i >= last_wvl:
                    if len(coeffs_per_wvl) > 0:
                        fit_prev = coeffs_per_wvl[-1]['fit']
                        I_full = data[i, 0, :, :]
                        for ch, idx in (('Q', 1), ('U', 2), ('V', 3)):
                            if ch not in channels or ch not in fit_prev:
                                continue
                            Y_full = data[i, idx, :, :]
                            Yhat_full = fit_prev[ch]['predict'](I_full)
                            M = np.isfinite(Y_full) & np.isfinite(Yhat_full)
                            out = data_corrected[i, idx, :, :]
                            out[M] = Y_full[M] - Yhat_full[M]
                    continue

                sub = data[i, :, ysl, xsl]  # (4, h, w)
                if dualI is not None:
                    I0 = dualI[0, i, ysl, xsl].astype(float)
                    I1 = dualI[1, i, ysl, xsl].astype(float)
                    dual = np.array([I0,I1])
                else:
                    dual = None

                (coefs, fit) = evaluate_crosstalk(sub, dualI = dual, **ec_kwargs)

                I_full = data[i, 0, :, :]
                for ch, idx in (('Q', 1), ('U', 2), ('V', 3)):
                    if ch not in channels or ch not in fit:
                        continue
                    Y_full = data[i, idx, :, :]
                    Yhat_full = fit[ch]['predict'](I_full)
                    M = np.isfinite(Y_full) & np.isfinite(Yhat_full)
                    out = data_corrected[i, idx, :, :] 
                    out[M] = Y_full[M] - Yhat_full[M]

                coeffs_per_wvl.append({
                    'wavelength_index': i,
                    'fit': fit,               
                    'tuple': coefs            
                })

                pixel_counts.append(int(np.isfinite(sub[0]).sum()))

            info.update({
                'mode': 'per_wavelength',
                'coeffs_global': None,
                'coeffs_per_wvl': coeffs_per_wvl,
                'pixel_counts': pixel_counts
            })
            return data_corrected, info

        # ======================== MÉTODO JAEGGLI / ALL_WVLS / CORR =========================
        elif method.lower() in ('jaeggli', 'all_wvls', 'corr'):
            # Optimizamos reordenamiento de ejes de moveaxis a transpose
            dat = np.transpose(data, (2, 3, 0, 1))  # (x, y, wl, st)

            Nwaves = dat.shape[2]
            full_data = dat.copy()
            # recorte espacial
            dat = dat[reg[0]:reg[1], reg[2]:reg[3], :, :]

            # recorte espectral para métrica
            if last_wvl is not None:
                dat = dat[:, :, :int(last_wvl), :]

            V_mean_corr = dat[:, :, :, 3] - np.mean(dat[:, :, 0, 3], axis=(0, 1))
            pmap = np.max(np.abs(V_mean_corr) / (dat[:, :, :, 0] + 1e-12), axis=2)
            notpolar = np.argwhere(pmap < pthresh)
            nyidx = notpolar[:, 0]
            nzidx = notpolar[:, 1]
            weak_region = dat[nyidx, nzidx, :, :] 

            if plots:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(6, 6))
                im = ax.imshow(pmap, cmap='viridis')
                ax.set_title('Fractional polarization map')
                plt.colorbar(im, ax=ax)
                plt.show(); plt.close()

            # ------------- Ajuste Jaeggli -------------
            if method.lower() == 'jaeggli':
                if MM1a is not None:
                    iMM1a = np.linalg.inv(MM1a)
                    data_inv = np.einsum('ij,abcj->abci', iMM1a, full_data)
                    data_corrected = np.transpose(data_inv, (2, 3, 0, 1))
                    info.update({'MM1a': MM1a})
                    return data_corrected, info
                else:
                    D0, theta0, chi0 = 0.5, 0.0, 0.0
                    res = minimize(lambda x: fitfunc1(x, weak_region),
                                x0=(D0, theta0, chi0),
                                options={'maxiter': 1000, 'disp': verbose})
                    MM1a = polmodel1(res.x[0], res.x[1], res.x[2])
                    iMM1a = np.linalg.inv(MM1a)
                    data_inv = np.einsum('ij,abcj->abci', iMM1a, full_data)
                    data_corrected = np.transpose(data_inv, (2, 3, 0, 1))
                    info.update({'MM1a': MM1a, 'optimizer_success': bool(res.success)})
                    return data_corrected, info

            # ------------- Ajuste por todas las λ -------------
            elif method.lower() in ('all_wvls', 'corr'):
                MM1a_arr = np.zeros((Nwaves, 4, 4), dtype=float)
                data_corr_xywl = full_data.copy()
                D0, theta0, chi0 = 0.5, 0.0, 0.0
                for wvli in range(Nwaves):
                    if verbose:
                        logging.info(f"[Jaeggli] Fitting wave {wvli+1}/{Nwaves}")
                    fun = lambda x: fitfunc1(x, weak_region, wvl=wvli, method=method)
                    res = minimize(fun, x0=(D0, theta0, chi0),
                                options={'maxiter': 1000, 'disp': verbose})
                    MM1 = polmodel1(res.x[0], res.x[1], res.x[2])
                    MM1a_arr[wvli, :, :] = MM1
                    iMM1 = np.linalg.inv(MM1)
                    data_inv = np.einsum('ij,abcj->abci', iMM1, full_data)
                    data_corr_xywl[:, :, wvli, :] = data_inv[:, :, wvli, :]

                data_corrected = np.transpose(data_corr_xywl, (2, 3, 0, 1))
                info.update({'MM1a': MM1a_arr})
                return data_corrected, info

        else:
            raise ValueError(f"Unknown method='{method}'. Use 'standard' | 'jaeggli' | 'all_wvls' | 'corr'.")

    elif strategy == 'sequential':
        # --- PASO 1: Ajuste principal ---
        # Se ejecutará 'jaeggli' o 'standard' dependiendo de lo que elijas en el parámetro 'method'
        logging.info(f"Iniciando PASO 1 secuencial: Corrección principal usando método '{method}'...")
        
        data_step1, info_step1 = fit_mueller_matrix(
            data, strategy='simultaneous', pthresh=pthresh, norm=norm,
            region=region, method=method, last_wvl=last_wvl, plots=plots,
            roi=roi, norm_wave=norm_wave, verbose=verbose, ctmethod=ctmethod,
            MM1a=MM1a, dualI=dualI, 
            use_local=False,           # Forzamos SIN gradientes en este primer paso
            # fit_grad=False,            # Ajuste completo
            local_order=local_order,
            deriv_sigma=deriv_sigma, local_ridge_lambda=local_ridge_lambda,
            aggregate_wavelengths=aggregate_wavelengths, channels=channels
        )

        # --- PASO 2: Ajuste de gradientes espaciales sobre el RESIDUO ---
        # Independientemente del método del Paso 1, el Paso 2 SIEMPRE usa 'standard' con derivadas
        logging.info("Iniciando PASO 2 secuencial: Corrección de gradientes ('standard' con derivadas)...")
        
        data_corrected, info_step2 = fit_mueller_matrix(
            data_step1, strategy='simultaneous', pthresh=pthresh, norm=False, # No normalizar el residuo
            region=region, method='standard', last_wvl=last_wvl, plots=False,
            roi=roi, norm_wave=norm_wave, verbose=verbose, ctmethod=ctmethod,
            MM1a=None, dualI=dualI, 
            use_local=True,            # Forzamos USO de gradientes
            # fit_grad=True,             # Forzamos a=0, b=0 para ajustar solo el residuo
            local_order=local_order,
            deriv_sigma=deriv_sigma, local_ridge_lambda=local_ridge_lambda,
            aggregate_wavelengths=aggregate_wavelengths, channels=channels
        )

        # Actualizamos la información de salida para saber exactamente qué se ha ejecutado
        info.update({
            'mode': 'sequential',
            'method_step1': method,
            'method_step2': 'standard_gradients',
            'step1_info': info_step1,
            'step2_info': info_step2
        })
        
        return data_corrected, info

def apply_crosstalk_coeffs_standard(
        data,
        coeffs,
        *,
        per_wavelength=False,        
        use_local=False,
        local_order=1,
        deriv_sigma=0.0,
        channels=('Q','U','V')
    ):
    """
    Aplica corrección de crosstalk I->(Q,U,V) usando COEFICIENTES DADOS por el usuario,
    sin calcular nada. SOLO modo STANDARD (lineal o local).
    """

    if data.ndim != 4 or data.shape[1] != 4:
        raise ValueError("data debe ser (wvl,4,x,y) con [I,Q,U,V]")

    n_wvl, _, X, Y = data.shape
    data_corr = np.copy(data)

    # -------------------------------------------------------------
    # Derivadas por λ (si use_local=True) -> PRECOMPUTACIÓN
    # -------------------------------------------------------------
    if use_local:
        Ix = np.zeros((n_wvl, X, Y), float)
        Iy = np.zeros((n_wvl, X, Y), float)
        L  = np.zeros((n_wvl, X, Y), float) if local_order >= 2 else None

        for i in range(n_wvl):
            Ib = data[i,0]
            if deriv_sigma > 0:
                Ib = gaussian_filter(Ib, deriv_sigma)
            Ix[i] = sobel(Ib, axis=1)
            Iy[i] = sobel(Ib, axis=0)
            if local_order >= 2:
                L[i] = laplace(Ib)

    # -------------------------------------------------------------
    # BUCLE PRINCIPAL POR λ
    # -------------------------------------------------------------
    for i in range(n_wvl):
        I = data[i,0]

        if per_wavelength:
            coeff_i = coeffs[i]
        else:
            coeff_i = coeffs

        for ch, idx in (('Q',1), ('U',2), ('V',3)):
            if ch not in channels:
                continue
            if ch not in coeff_i:
                continue

            a = float(coeff_i[ch].get('a',0.0))
            b = float(coeff_i[ch].get('b',0.0))
            c = float(coeff_i[ch].get('c',0.0))
            d = float(coeff_i[ch].get('d',0.0))
            e = float(coeff_i[ch].get('e',0.0))

            if not use_local:
                Yhat = a + b * I
            else:
                Yhat = a + b*I + c*Ix[i] + d*Iy[i]
                if local_order >= 2:
                    Yhat = Yhat + e * L[i]

            Y = data[i, idx]
            M = np.isfinite(Y) & np.isfinite(Yhat)
            data_corr[i, idx][M] = Y[M] - Yhat[M]

    info = dict(
        applied=True,
        per_wavelength=bool(per_wavelength),
        use_local=bool(use_local),
        local_order=int(local_order),
        deriv_sigma=float(deriv_sigma),
        channels=channels
    )
    return data_corr, info

def fit_mueller_matrix_tiled(
    data: np.ndarray,
    strategy: str = "simultaneous",
    pthresh: float = 0.02,
    norm: bool = False,
    quadrants: int = 8,
    overlap_fraction: float = 0.25,          # NUEVO: % de solapamiento para borrar líneas (0.25 = 25%)
    region: List[int] = [0, -1, 0, -1],      
    method: str = 'standard',            
    last_wvl: Optional[int] = None,
    plots: bool = False,
    roi: List[int] = [0, -1, 0, -1],         
    norm_wave: int = -1,
    verbose: bool = False,
    ctmethod: str = 'linfit',
    MM1a: Optional[np.ndarray] = None,
    # ===== NUEVO: Mismas capacidades que la función base =====
    dualI: Optional[np.ndarray] = None,
    use_local: bool = False,
    # fit_grad: bool = False,
    local_order: int = 1,
    deriv_sigma: float = 0.0,
    local_ridge_lambda: float = 0.0,
    aggregate_wavelengths: bool = False,
    channels: Tuple[str, ...] = ('Q','U','V')
) -> Tuple[np.ndarray, List[Tuple[int,int,int,int]], List[Dict[str,Any]]]:
    """
    Divide la imagen en una rejilla con solapamiento (feathering) para evitar líneas 
    de costura, y aplica fit_mueller_matrix en cada zona.
    """
    
    if data.ndim != 4 or data.shape[1] != 4:
        raise ValueError("Expected data shape (wavelength, 4, x, y) with Stokes=[I,Q,U,V].")

    H, W = data.shape[-2], data.shape[-1]
    data_float = data.astype(float, copy=False)

    # --- 1. Normalización GLOBAL (Para evitar saltos de brillo entre tiles) ---
    if norm:
        y1_roi, y2_roi, x1_roi, x2_roi = roi
        y2_roi = H if y2_roi < 0 else y2_roi
        x2_roi = W if x2_roi < 0 else x2_roi
        I0 = data_float[norm_wave, 0, y1_roi:y2_roi, x1_roi:x2_roi]
        nf = float(np.nanmedian(I0))
        if nf > 0 and np.isfinite(nf):
            data_float = data_float / nf
    
    # --- 2. Generador de cuadrantes con solapamiento ---
    def _generate_tiles_with_overlap(H, W, n, overlap_frac):
        base_h = H / n
        base_w = W / n
        pad_h = int(round(base_h * overlap_frac))
        pad_w = int(round(base_w * overlap_frac))
        
        tiles = []
        for i in range(n):
            for j in range(n):
                # Expandimos el tile por el pad, sin salirnos de la imagen
                y1 = max(0, int(round(i * base_h)) - pad_h)
                y2 = min(H, int(round((i + 1) * base_h)) + pad_h)
                x1 = max(0, int(round(j * base_w)) - pad_w)
                x2 = min(W, int(round((j + 1) * base_w)) + pad_w)
                
                if y2 > y1 and x2 > x1:
                    tiles.append((y1, y2, x1, x2, pad_h, pad_w))
        return tiles

    # --- 3. Ventana de atenuación lineal (Bilinear Tapering) ---
    def _make_taper_window(y1, y2, x1, x2, H, W, pad_h, pad_w):
        """Crea una ventana 2D que cae a 0 en los bordes solapados para un fundido perfecto."""
        h, w = y2 - y1, x2 - x1
        wy = np.ones(h, dtype=float)
        wx = np.ones(w, dtype=float)
        
        # Atenuar borde superior (solo si no es el borde absoluto de la imagen)
        if y1 > 0 and pad_h > 0:
            wy[:pad_h] = np.linspace(0, 1, pad_h)
        # Atenuar borde inferior
        if y2 < H and pad_h > 0:
            wy[-pad_h:] = np.linspace(1, 0, pad_h)
        # Atenuar borde izquierdo
        if x1 > 0 and pad_w > 0:
            wx[:pad_w] = np.linspace(0, 1, pad_w)
        # Atenuar borde derecho
        if x2 < W and pad_w > 0:
            wx[-pad_w:] = np.linspace(1, 0, pad_w)
            
        return wy[:, None] * wx[None, :]

    tiles_info = _generate_tiles_with_overlap(H, W, quadrants, overlap_fraction)
    
    # Acumuladores para la reconstrucción suave
    data_corr_accum = np.zeros_like(data_float)
    weight_accum = np.zeros((H, W), dtype=float)
    infos: List[Dict[str, Any]] = []

    # --- 4. Bucle principal de procesado por Tile ---
    for q, (y1, y2, x1, x2, pad_h, pad_w) in enumerate(tiles_info):
        sub_data = data_float[:, :, y1:y2, x1:x2]
        
        # Trocear también el haz dual si existe
        if dualI is not None:
            sub_dualI = dualI[..., y1:y2, x1:x2]
        else:
            sub_dualI = None

        # Llamada al núcleo matemático (Forzamos norm=False porque ya se hizo globalmente)
        data_corr_tile, info = fit_mueller_matrix(
            sub_data,
            strategy=strategy,
            pthresh=pthresh,
            norm=False,                
            region=[0, -1, 0, -1],     
            method=method,
            last_wvl=last_wvl,
            plots=plots,
            verbose=verbose,
            ctmethod=ctmethod,
            MM1a=MM1a,
            dualI=sub_dualI,
            use_local=use_local,
            # fit_grad=fit_grad,
            local_order=local_order,
            deriv_sigma=deriv_sigma,
            local_ridge_lambda=local_ridge_lambda,
            aggregate_wavelengths=aggregate_wavelengths,
            channels=channels
        )

        # Crear máscara de fundido y acumular
        window = _make_taper_window(y1, y2, x1, x2, H, W, pad_h, pad_w)
        
        data_corr_accum[:, :, y1:y2, x1:x2] += data_corr_tile * window[None, None, :, :]
        weight_accum[y1:y2, x1:x2] += window

        info = dict(info)  
        info['tile_bounds'] = (y1, y2, x1, x2)
        infos.append(info)

        if verbose:
            logging.info(f"[Tile {q+1}/{len(tiles_info)}] OK -> bounds=({y1}:{y2}, {x1}:{x2})")

    # --- 5. Promediado final ponderado ---
    # Dividimos entre la suma de pesos para normalizar el fundido (evitando div/0)
    weight_accum = np.maximum(weight_accum, 1e-12)
    data_corr_final = data_corr_accum / weight_accum[None, None, :, :]

    return data_corr_final, [(t[0], t[1], t[2], t[3]) for t in tiles_info], infos

def ver_planos_ajuste(planes, samples):
    def _pts(pts):
        if pts is None:
            return []
        if isinstance(pts, (list, tuple)):
            if len(pts) == 0: return []
            if isinstance(pts[0], (list, tuple, np.ndarray)) and len(pts[0]) == 3:
                return [tuple(map(float, p)) for p in pts]
            if len(pts) == 3 and all(isinstance(v, (int,float,np.floating)) for v in pts):
                return [tuple(map(float, pts))]
            return []
        if isinstance(pts, np.ndarray):
            pts = np.asarray(pts)
            if pts.ndim == 2 and pts.shape[1] == 3:
                return [tuple(map(float, p)) for p in pts]
            if pts.ndim == 1 and pts.size == 3:
                return [tuple(map(float, pts))]
        return []

    any_plot = False
    for ch in sorted(planes.keys()):
        for t in sorted(planes[ch].keys()):

            coefs = planes[ch][t]
            if len(coefs) == 3:
                A, B, C = coefs
                Dxy = 0.0
            else:
                A, B, C, Dxy = coefs

            pts = _pts(samples.get(ch, {}).get(t, []))

            print(f"\n=== {ch} / {t} ===")
            print(f"Plano: z = {A:.6e} * x + {B:.6e} * y + {C:.6e}")
            if not pts:
                print("Puntos: (ninguno)")
                continue

            print("Puntos (cx, cy, z):")
            for (x,y,z) in pts:
                print(f"  ({x:.1f}, {y:.1f}, {z:.6e})")

            if len(pts) < 3:
                print(f"Aviso: solo {len(pts)} punto(s). No se dibuja el plano.")
                continue

            xs = np.array([p[0] for p in pts], float)
            ys = np.array([p[1] for p in pts], float)
            zs = np.array([p[2] for p in pts], float)

            xx = np.linspace(xs.min(), xs.max(), 20)
            yy = np.linspace(ys.min(), ys.max(), 20)
            XX, YY = np.meshgrid(xx, yy)
            ZZ = A*XX + B*YY + C + Dxy*(XX*YY)

            fig = plt.figure(figsize=(5,4))
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(xs, ys, zs, c='r', s=40, label='Puntos')
            ax.plot_surface(XX, YY, ZZ, alpha=0.5, cmap='viridis')
            ax.set_title(f"{ch} – {t}")
            ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("coef")
            ax.legend()
            plt.tight_layout(); plt.show()
            any_plot = True

    if not any_plot:
        print("\nNo se ha dibujado ninguna superficie (faltan términos con >=3 puntos).")

def fit_mueller_matrix_rois_with_plane(data, rois, channels=('Q','U','V'), **kw):
    L, _, Y, X = data.shape
    name2idx = {'I':0,'Q':1,'U':2,'V':3}
    TERMS = ('a','b','c','d','e')

    def recorte(cx,cy,lado):
        h = lado//2
        x1=max(0,cx-h); x2=min(X,cx+h+(lado%2))
        y1=max(0,cy-h); y2=min(Y,cy+h+(lado%2))
        return x1,x2,y1,y2

    def extraer_coef(info):
        out={ch:{} for ch in channels}
        per = info.get('coeffs_per_wvl',[])
        if not per: return {}
        for ch in channels:
            acc={t:[] for t in TERMS}
            for elt in per:
                fd = elt['fit'].get(ch,{})
                for t in TERMS:
                    if t in fd: acc[t].append(float(fd[t]))
            out[ch] = {t: np.nanmean(acc[t]) for t in TERMS if acc[t]}
        return out

    samples = {ch:{t:[] for t in TERMS} for ch in channels}
    infos=[]

    for (cx,cy,lado) in rois:
        print(f"Procesando ROI centrada en ({cx},{cy}) con lado {lado}...")
        x1,x2,y1,y2 = recorte(cx,cy,lado)
        sub = data[:,:,y1:y2, x1:x2]

        _, info = fit_mueller_matrix(sub, channels=channels, **kw)
        infos.append(info)

        coefs = extraer_coef(info)
        for ch in channels:
            for t in TERMS:
                if t in coefs.get(ch,{}):
                    samples[ch][t].append((cx,cy,float(coefs[ch][t])))
        logging.info(f"ROI ({cx},{cy}) -> coeficientes extraídos: {coefs}")

    planes={ch:{} for ch in channels}
    maps  ={ch:{} for ch in channels}

    xs=np.arange(X)[None,:]
    ys=np.arange(Y)[:,None]

    for ch in channels:
        for t in TERMS:
            pts=samples[ch][t]
            if len(pts)<3: continue

            P = np.array([[x, y, x*y, 1.0] for (x,y,_) in pts], float)
            z = np.array([v for (_,_,v) in pts], float)
            A, B, Dxy, C = np.linalg.lstsq(P, z, rcond=None)[0]
            planes[ch][t] = (A, B, C, Dxy)             
            maps[ch][t]   = A*xs + B*ys + C + Dxy*(ys*xs)  

    data_corr = data.copy()
    for l in range(L):
        I = data_corr[l,0]
        dIx = sobel(I, axis=1)
        dIy = sobel(I, axis=0)

        for ch in channels:
            j = name2idx[ch]
            CT = np.zeros_like(I)
            for t,(arr) in maps[ch].items():
                if t=='a': CT += arr
                if t=='b': CT += arr * I
                if t=='c': CT += arr * dIx
                if t=='d': CT += arr * dIy
            data_corr[l,j] -= CT

    return data_corr, planes, maps, infos, samples

def write_crosstalk_header(
    header: "Header",
    info: dict,
    channels=('Q','U','V'),
    index_width: int = 3,   
):
    header.set('CROSTALK', 1, 'Was crosstalk correction applied?')
    header.set('CROSMETH', str(info.get('method','standard')), 'Crosstalk method')
    header.set('CROSMODE', str(info.get('mode','per_wavelength')), 'Fit mode: per_wavelength or global')

    local_flag = 0
    local_order = 0
    deriv_sigma = 0.0
    ridge_lambda = 0.0

    if info.get('mode') == 'global':
        fitG = info.get('coeffs_global', {})
        for ch in channels:
            if ch in fitG:
                local_flag = 1 if fitG[ch].get('model','linear').startswith('local') else 0
                if local_flag:
                    try:
                        local_order = int(fitG[ch].get('model','linear')[-1])
                    except Exception:
                        local_order = 1
                deriv_sigma = float(fitG[ch].get('deriv_sigma', 0.0))
                ridge_lambda = float(fitG[ch].get('lambda', 0.0))
                break
    else:
        per = info.get('coeffs_per_wvl', [])
        if len(per) > 0:
            fit0 = per[0].get('fit', {})
            for ch in channels:
                if ch in fit0:
                    local_flag = 1 if fit0[ch].get('model','linear').startswith('local') else 0
                    if local_flag:
                        try:
                            local_order = int(fit0[ch].get('model','linear')[-1])
                        except Exception:
                            local_order = 1
                    deriv_sigma = float(fit0[ch].get('deriv_sigma', 0.0))
                    ridge_lambda = float(fit0[ch].get('lambda', 0.0))
                    break

    header.set('CROSMLOC', int(local_flag), 'Local model (derivatives) used?')
    header.set('CROSLORD', int(local_order), 'Local model order (1 or 2)')
    header.set('CROSDSIG', float(deriv_sigma), 'Derivative pre-smoothing sigma')
    header.set('CROSRIDG', float(ridge_lambda), 'Ridge lambda (local model)')

    mode = info.get('mode','per_wavelength')

    def _write_coeffs_for_channel(prefix: str, a: float, b: float, c: float, d: float, e: float, suffix: str):
        header.set(f'S{prefix}_{suffix}', float(b), f'slope b for {prefix}')
        header.set(f'C{prefix}_{suffix}', float(a), f'intercept a for {prefix}')
        header.set(f'DX{prefix}_{suffix}', float(c), f'd/dx coeff c for {prefix}')
        header.set(f'DY{prefix}_{suffix}', float(d), f'd/dy coeff d for {prefix}')
        header.set(f'DL{prefix}_{suffix}', float(e), f'laplacian coeff e for {prefix}')

    if mode == 'global':
        fitG = info.get('coeffs_global', {})
        for ch in channels:
            fd = fitG.get(ch, {})
            a = float(fd.get('a', 0.0))
            b = float(fd.get('b', 0.0))
            c = float(fd.get('c', 0.0))
            d = float(fd.get('d', 0.0))
            e = float(fd.get('e', 0.0))
            _write_coeffs_for_channel(ch, a, b, c, d, e, 'G')
    else:
        per = info.get('coeffs_per_wvl', [])
        if not per:
            return header
        for dct in per:
            widx = int(dct.get('wavelength_index', -1))
            suffix = str(widx).zfill(index_width) if widx >= 0 else 'UNK'
            fitW = dct.get('fit', {})
            for ch in channels:
                fd = fitW.get(ch, {})
                a = float(fd.get('a', 0.0))
                b = float(fd.get('b', 0.0))
                c = float(fd.get('c', 0.0))
                d = float(fd.get('d', 0.0))
                e = float(fd.get('e', 0.0))
                _write_coeffs_for_channel(ch, a, b, c, d, e, suffix)

    return header

def _grid_tiles(H: int, W: int, y1: int, y2: int, x1: int, x2: int, divisions: int) -> List[Tuple[int,int,int,int]]:
    h = max(0, y2 - y1); w = max(0, x2 - x1)
    if h <= 0 or w <= 0:
        raise ValueError("Region must have positive area.")
    dh = h / divisions; dw = w / divisions
    tiles = []
    for i in range(divisions):
        for j in range(divisions):
            yt1 = int(round(y1 + i * dh)); yt2 = int(round(y1 + (i+1) * dh))
            xt1 = int(round(x1 + j * dw)); xt2 = int(round(x1 + (j+1) * dw))
            yt1 = max(0, min(H, yt1)); yt2 = max(0, min(H, yt2))
            xt1 = max(0, min(W, xt1)); xt2 = max(0, min(W, xt2))
            if yt2 > yt1 and xt2 > xt1:
                tiles.append((yt1, yt2, xt1, xt2))
    return tiles

def _pearson(u, v):
    u = np.asarray(u).ravel(); v = np.asarray(v).ravel()
    m = np.isfinite(u) & np.isfinite(v)
    if m.sum() < 3:
        return np.nan
    u = u[m]; v = v[m]
    su = np.std(u); sv = np.std(v)
    if su == 0 or sv == 0:
        return np.nan
    return float(np.corrcoef(u, v)[0,1])

def _linfit_clip(x, y, n_sigma=3):
    x = x.ravel(); y = y.ravel()
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan, np.nan
    xv = x[m]; yv = y[m]
    b1, a1 = np.polyfit(xv, yv, 1)
    res = yv - (a1 + b1*xv)
    s = np.std(res) + 1e-12
    keep = np.abs(res) < n_sigma * s
    if keep.sum() < 3:
        return float(a1), float(b1)
    b2, a2 = np.polyfit(xv[keep], yv[keep], 1)
    return float(a2), float(b2)

def _ridge_solve(A, y, lam=0.0):
    AtA = A.T @ A
    Aty = A.T @ y
    if lam > 0:
        L = np.eye(AtA.shape[0])
        L[0,0] = 0.0
        AtA = AtA + lam * L
    return np.linalg.solve(AtA, Aty)

def fit_interference_Iref_to_Q_tiled(
    data: np.ndarray,                      
    ref_wvl: Optional[int] = -1,           
    divisions: int = 8,
    region: List[int] = [0,-1,0,-1],       
    ctmethod: str = 'linfit',              
    n_sigma: float = 3,                    
    pthresh_intensity: float = 0.0,        
    use_local: bool = False,               
    local_order: int = 1,                  
    deriv_sigma: float = 0.0,              
    local_ridge_lambda: float = 0.0,       
    apply: bool = True,                    
    show_grid: bool = False,
    verbose: bool = False
) -> Tuple[np.ndarray, List[Tuple[int,int,int,int]], Dict[str, Any]]:
    
    if data.ndim != 4 or data.shape[1] != 4:
        raise ValueError("Expected data shape (wavelength, 4, x, y) with Stokes=[I,Q,U,V].")

    n_wvl, _, H, W = data.shape
    if ref_wvl is None:
        ref_wvl = -1
    if ref_wvl < 0:
        ref_wvl = n_wvl + ref_wvl  
    ref_wvl = int(np.clip(ref_wvl, 0, n_wvl-1))

    y1, y2, x1, x2 = region
    if y2 < 0: y2 = H
    if x2 < 0: x2 = W
    tiles = _grid_tiles(H, W, y1, y2, x1, x2, divisions)

    Iref = data[ref_wvl, 0, :, :].astype(float)
    Iref_b = gaussian_filter(Iref, sigma=deriv_sigma) if (use_local and deriv_sigma > 0) else Iref
    Ix_ref = sobel(Iref_b, axis=1) if use_local else None
    Iy_ref = sobel(Iref_b, axis=0) if use_local else None
    L_ref  = laplace(Iref_b)       if (use_local and local_order >= 2) else None

    n_tiles = len(tiles)
    a = np.full((n_wvl, n_tiles), np.nan, float)
    b = np.full((n_wvl, n_tiles), np.nan, float)
    c = np.full((n_wvl, n_tiles), 0.0, float)
    d = np.full((n_wvl, n_tiles), 0.0, float)
    e = np.full((n_wvl, n_tiles), 0.0, float)
    corr = np.full((n_wvl, n_tiles), np.nan, float)

    data_corr = np.copy(data)

    if show_grid:
        I0 = data[0, 0]
        fig, ax = plt.subplots(figsize=(7,7))
        im = ax.imshow(I0, cmap='gray')
        for (yt1, yt2, xt1, xt2) in tiles:
            ax.add_patch(Rectangle((xt1, yt1), xt2-xt1, yt2-yt1, fill=False, edgecolor='r', lw=0.8))
        ax.set_title(f"Tiles {divisions}×{divisions} (I @ λref={ref_wvl})")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        plt.tight_layout(); plt.show()

    for t_idx, (yt1, yt2, xt1, xt2) in enumerate(tiles):
        X  = Iref[yt1:yt2, xt1:xt2]
        if use_local:
            Xx = Ix_ref[yt1:yt2, xt1:xt2]
            Xy = Iy_ref[yt1:yt2, xt1:xt2]
            Lx = L_ref[yt1:yt2, xt1:xt2] if local_order >= 2 else None

        M0 = np.isfinite(X)
        if pthresh_intensity > 0:
            M0 &= (X > pthresh_intensity)

        for w in range(n_wvl):
            Y = data[w, 1, yt1:yt2, xt1:xt2].astype(float)  
            M = M0 & np.isfinite(Y)
            if M.sum() < 5:
                continue

            corr[w, t_idx] = _pearson(X[M], Y[M])

            if not use_local:
                if ctmethod == 'linfit':
                    a_w, b_w = _linfit_clip(X[M], Y[M], n_sigma=n_sigma)
                elif ctmethod == 'rms':
                    xv = X[M].ravel(); yv = Y[M].ravel()
                    denom = (xv**2).sum()
                    b_w = float((xv @ yv) / denom) if denom > 0 else np.nan
                    a_w = float(yv.mean() - b_w * xv.mean()) if np.isfinite(b_w) else np.nan
                else:
                    raise ValueError("ctmethod must be 'linfit' or 'rms'.")
                a[w, t_idx], b[w, t_idx] = a_w, b_w
                if apply and np.isfinite(a_w) and np.isfinite(b_w):
                    Yhat = a_w + b_w * X
                    Mf = np.isfinite(Yhat) & np.isfinite(Y)
                    data_corr[w, 1, yt1:yt2, xt1:xt2][Mf] = Y[Mf] - Yhat[Mf]

            else:
                cols = [np.ones_like(X), X, Xx, Xy]
                if local_order >= 2:
                    cols.append(Lx)
                mats = [col[M].ravel() for col in cols]
                A = np.column_stack(mats)
                yv = Y[M].ravel()
                theta = _ridge_solve(A, yv, lam=float(local_ridge_lambda))
                a_w, b_w, c_w, d_w = theta[:4]
                e_w = float(theta[4]) if (local_order >= 2 and theta.size >= 5) else 0.0
                a[w, t_idx], b[w, t_idx], c[w, t_idx], d[w, t_idx], e[w, t_idx] = \
                    float(a_w), float(b_w), float(c_w), float(d_w), float(e_w)

                if apply and np.all(np.isfinite([a_w, b_w, c_w, d_w, e_w])):
                    Yhat = a_w + b_w*X + c_w*Xx + d_w*Xy + (e_w*Lx if local_order >= 2 else 0.0)
                    Mf = np.isfinite(Yhat) & np.isfinite(Y)
                    data_corr[w, 1, yt1:yt2, xt1:xt2][Mf] = Y[Mf] - Yhat[Mf]

        if verbose:
            print(f"[tile {t_idx+1}/{n_tiles}] done.")

    out = {
        'ref_wvl': ref_wvl,
        'divisions': divisions,
        'region': (y1,y2,x1,x2),
        'use_local': bool(use_local),
        'local_order': int(local_order),
        'deriv_sigma': float(deriv_sigma),
        'ridge': float(local_ridge_lambda),
        'a': a, 'b': b, 'c': c, 'd': d, 'e': e,
        'corr': corr,                       
        'tiles': tiles
    }
    return data_corr, tiles, out

def write_interference_Iref_to_Q_header(
    header: "Header",
    a: np.ndarray,        
    b: np.ndarray,        
    c: np.ndarray = None, 
    d: np.ndarray = None, 
    e: np.ndarray = None, 
    *,
    ref_wvl: int,
    divisions: int,
    tiles: list,                   
    use_local: bool = False,
    local_order: int = 1,
    deriv_sigma: float = 0.0,
    ridge_lambda: float = 0.0,
    index_width: int = 3,          
    after_keyword: str = 'ALIGMETH'
):
    header.set('CROSTALK', 1, 'Was crosstalk correction applied?', after=after_keyword if after_keyword in header else None)
    header.set('XK_MODE', 'IREF2Q', 'Interference mode: I(ref) -> Q')
    header.set('XK_REF', int(ref_wvl), 'Reference wavelength index used for I')
    header.set('XK_DIV', int(divisions), 'Grid divisions per axis (n x n)')
    header.set('XK_NTIL', int(len(tiles)), 'Number of tiles')
    header.set('XK_LCL', int(bool(use_local)), 'Local model (derivatives) used?')
    header.set('XK_ORD', int(local_order), 'Local model order (1 or 2)')
    header.set('XK_DSIG', float(deriv_sigma), 'Derivative pre-smoothing sigma')
    header.set('XK_RIDG', float(ridge_lambda), 'Ridge lambda (local model)')

    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n_wvl, n_tiles = a.shape

    def _nz(x, shape):
        if x is None:
            return np.zeros(shape, dtype=float)
        xx = np.asarray(x, dtype=float)
        if xx.shape != shape:
            raise ValueError(f"Shape mismatch: expected {shape}, got {xx.shape}")
        return xx

    c = _nz(c, a.shape)
    d = _nz(d, a.shape)
    e = _nz(e, a.shape)

    for iw in range(n_wvl):
        for it in range(n_tiles):
            suf = f"{iw:0{index_width}d}{it:0{index_width}d}"   
            header[f'IA{suf}'] = (float(a[iw, it]), f'Interference I->Q: intercept a (w={iw},t={it})')
            header[f'IB{suf}'] = (float(b[iw, it]), f'Interference I->Q: slope b (w={iw},t={it})')
            header[f'IC{suf}'] = (float(c[iw, it]), f'Interference I->Q: d/dx c (w={iw},t={it})')
            header[f'ID{suf}'] = (float(d[iw, it]), f'Interference I->Q: d/dy d (w={iw},t={it})')
            header[f'IL{suf}'] = (float(e[iw, it]), f'Interference I->Q: laplacian e (w={iw},t={it})')

    return header

def load_polynomial_coeffs_csv(csv_path):
    with open(csv_path) as f:
        r = csv.DictReader(f)
        rows = list(r)

    coeffs = {}
    for row in rows:
        iw = int(row["wavelength"])
        degx = int(row["degx"])
        degy = int(row["degy"])

        surfA = {}
        surfB = {}

        for key,val in row.items():
            if key.startswith("a_"):
                _,i,j = key.split("_")
                surfA[(int(i),int(j))] = float(val)
            if key.startswith("b_"):
                _,i,j = key.split("_")
                surfB[(int(i),int(j))] = float(val)

        coeffs[iw] = {
            "degx": degx,
            "degy": degy,
            "a": surfA,
            "b": surfB
        }
    return coeffs

def normalize_pixel(x, y, nx, ny):
    xn = 2*(x+0.5)/nx - 1.0
    yn = 2*(y+0.5)/ny - 1.0
    return xn, yn

def eval_poly2d(x, y, surf):
    value = 0.0
    for (i,j), c in surf.items():
        value += c * (x**i) * (y**j)
    return value

def apply_polynomial_Iref2Q(Q, Iref, coeffs, wavelength):
    surfA = coeffs[wavelength]["a"]
    surfB = coeffs[wavelength]["b"]

    ny,nx = Q.shape
    Qcorr = np.zeros_like(Q)

    for y in range(ny):
        for x in range(nx):
            xn, yn = normalize_pixel(x, y, nx, ny)
            a = eval_poly2d(xn, yn, surfA)
            b = eval_poly2d(xn, yn, surfB)
            Qcorr[y,x] = Q[y,x] - (a + b * Iref[y,x])

    return Qcorr
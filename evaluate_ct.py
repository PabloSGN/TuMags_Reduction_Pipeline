import logging
import numpy as np
from scipy.ndimage import sobel, laplace, gaussian_filter
from scipy.stats import linregress

def evaluate_crosstalk(
        data,
        pthresh=0.05,
        region=[0,-1,0,-1],
        pthresh_intensity=0,
        channels=('Q','U','V'),
        use_local=False,
        local_order=1,
        deriv_sigma=0.0,
        local_ridge_lambda=0.0,
        fit_ab=True     
):
    """
    Crosstalk evaluation con opción de:
      - ajuste clásico completo (a,b,...)
      - ajuste local SOLO derivadas (a=0,b=1)
    """


    # ===============================
    # Helpers
    # ===============================

    def linfit_clip(x, y, n_sigma=5):
        slope, intercept, *_ = linregress(x, y)
        resid = y - (slope*x + intercept)
        std = np.std(resid) + 1e-12
        mask = np.abs(resid) < (n_sigma * std)
        slope2, intercept2, *_ = linregress(x[mask], y[mask])
        return intercept2, slope2

    def ridge_solve(A, y, lam=0.0):
        AtA = A.T @ A
        Aty = A.T @ y
        if lam > 0:
            L = np.eye(AtA.shape[0])
            L[0, 0] = 0.0
            AtA = AtA + lam * L
        return np.linalg.solve(AtA, Aty)

    def flatten_mask(mask, *arrays):
        idx = np.where(mask)
        return [arr[idx].ravel() for arr in arrays]

    def build_derivs(X):
        Xb = gaussian_filter(X, deriv_sigma) if deriv_sigma > 0 else X
        Ix = sobel(Xb, axis=1)
        Iy = sobel(Xb, axis=0)
        L = laplace(Xb) if local_order >= 2 else None
        return Ix, Iy, L

    # ===============================
    # Datos
    # ===============================

    y1,y2,x1,x2 = region
    sub = data[:, y1:y2, x1:x2].astype(float)

    I, Q, U, V = sub

    X2d = I   # simplificación: solo I en local

    comps = {'Q':Q, 'U':U, 'V':V}
    comps = {k:v for k,v in comps.items() if k in channels}

    eps = 1e-12
    V_mc = V - np.mean(V)

    sel = np.abs(V_mc)/(I+eps) < pthresh

    if pthresh_intensity != 0:
        sel &= I > pthresh_intensity

    results = {}

    # ===============================
    # LOOP canales
    # ===============================

    for name, Y in comps.items():

        if not use_local:
            # ----------------------------------
            # LINEAR NORMAL
            # ----------------------------------
            x, y = flatten_mask(sel, X2d, Y)

            if len(x) < 5:
                results[name]=(np.nan, np.nan,0,0,0)
                continue

            a,b = linfit_clip(x,y)

            results[name]=(a,b,0,0,0)
            continue

        # ======================================
        # LOCAL MODE
        # ======================================

        Ix, Iy, L = build_derivs(X2d)

        if fit_ab:
            # ----------------------------------
            # ORIGINAL LOCAL (fit completo)
            # ----------------------------------
            cols = [np.ones_like(X2d), X2d, Ix, Iy]
            if local_order >= 2:
                cols.append(L)

            mats = flatten_mask(sel, *cols, Y)
            *Xs, yv = mats

            A = np.column_stack([x.ravel() for x in Xs])
            theta = ridge_solve(A, yv, local_ridge_lambda)

            a,b,c,d = theta[:4]
            e = theta[4] if local_order>=2 else 0

        else:
            # ==================================
            # NO ajustar a y b
            # ==================================
            # Asumimos:
            #   a = 0
            #   b = 1
            # Datos YA corregidos
            #
            # Modelo:
            #   Y ≈ X + c·Ix + d·Iy
            # => residual = Y - X
            # ==================================

            residual = Y - X2d

            cols = [Ix, Iy]
            if local_order >= 2:
                cols.append(L)

            mats = flatten_mask(sel, *cols, residual)

            if len(mats) < 2:
                results[name]=(0.0,1.0,0.0,0.0,0.0)
                continue

            *Xs, yv = mats
            A = np.column_stack([x.ravel() for x in Xs])

            theta = ridge_solve(A, yv, local_ridge_lambda)

            c = theta[0]
            d = theta[1] if len(theta)>1 else 0
            e = theta[2] if local_order>=2 and len(theta)>2 else 0

            a = 0.0
            b = 1.0

        results[name]=(float(a),float(b),float(c),float(d),float(e))

    # ===============================
    # OUTPUT clásico
    # ===============================

    def unpack(k):
        return results.get(k,(np.nan,)*5)

    aQ,bQ,_,_,_ = unpack('Q')
    aU,bU,_,_,_ = unpack('U')
    aV,bV,_,_,_ = unpack('V')

    return aQ, aU, aV, bQ, bU, bV
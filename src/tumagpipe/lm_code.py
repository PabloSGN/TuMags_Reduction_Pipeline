from math import tau
import numpy as np
from matplotlib import pyplot as plt
import time
from copy import copy, deepcopy
from scipy.linalg import solve,LinAlgError

pardata = {'delta1': np.array([225, 225, 315, 315]),
            'delta2' : np.array([234, 125.26, 54.74, 305.26]),
            'theta1' : 0,
            'theta2' : 45,
            'pol_angle' : 0,
            'rot_inst' : 0}

"""
leastsq.py

author: Colin Clement
date: 2019-03-26

This module defines a Levenberg Marquardt optimization routine, featuring
geodesic acceleration from Mark Transtrum's PRE paper

    @article{transtrum2012geodesic,
      title={Geodesic acceleration and the small-curvature approximation
         for nonlinear least squares},
      author={Transtrum, Mark K and Sethna, James P},
      journal={arXiv preprint arXiv:1207.4999},
      year={2012}
    }

Usage:
    >>> opt = LM(residual_function, jacobian_function, func_args)
    >>> opt.leastsq(initial_point)

"""

class LM(object):
    def __init__(self, res, jac, args = ()):
        r"""
        Custom implementation of Levenburg-Marqardt. Solves the problem

            argmin_\theta \sum_i r_i(\theta)^2

        Parameters:
        -----------
            res: function
                residual function, takes parameters theta (arrays of length M)
                and optionally specified args, and returns residuals r_i
                of length N
            jac: function
                jacobian function which takes parameters theta (arrays of length
                M) and optionally specified args, and returns a jacobian matrix
                of shape (M, N)
            (optional)
            args: tuple of positional arguments to be passed into res
                and jac
        Usage:
        ------
        >>> opt = LM(res, jac)
        >>> opt.leastsq(initial_point)
        """
        self.res = res
        self.jac = jac
        self.args = args

        self.msgdict = {0: 'Failed number of steps reached maxtries',
                        1: 'Maximum number of iterations reached',
                        2: 'Convergence criterion satisfied'}

    def leastsq(self, p0, maxiter=40, delta=1E-6, accept=10.,
                reject=5., iprint=0, geo=True, **kwargs):
        """
        A custom implementation of the Levenburg Marquardt least-squares
        optimization algorithm.
        Parameters
        ----------
            p0: array_like
                initial list of parameter values
            (optional)
            maxiter: int
                maximum number of iterations
            delta: float
                finishing condition, algorithm quits if logprobability
                changes by less than delta
            accept: float
                factor to increase lambda upon accepting a step
            reject: float
                factor to decrease lambda upon rejecting a step
            geo: True/False
                whether or not to use geodesic acceleration
            iprint: integer
                0 for no messages, 1 for few messages, 2 for more messages
            kwargs:
                lamb : float
                    initial regularization parameter value
                maxtries: int
                    number of failed steps after which the
                    algorithm gives up reducing lambda
                solver: str
                    either 'cg' for conjugate gradient
                    or 'spsolver' for direct solver (faster but uses
                    ENORMOUS amounts of memory for standard scipy spsolve)
        returns:
            out_p, jtj, message : tuple
                opt_p (array of shape N_params), JTJ (Hessian), message (int)
                    message = 0 if number of failed steps is maxtries
                    message = 1 if the algorithm reach max iterations
                    message = 2 is algorithm achieved convergence criterion
        """
        lamb = float(kwargs.get('lamb', 10.)) #start with small downward grad steps
        h = kwargs.get('h', 1E-4)
        maxtries = kwargs.get('maxtries', 200)
        truncerr = 0.

        p1 = copy(p0)
        for itn in range(maxiter):
            r0 = self.res(p1, *self.args)
            nlnprob0 = r0.dot(r0)/2.

            if iprint:
                print("Itn {}: nlnprob = {}".format(itn, nlnprob0))

            j = self.jac(p1, *self.args)
            jtj = j.T.dot(j)
            jtr = j.T.dot(r0)
            
            success, tries, message = False, 0, 0
            while not success:
                if tries == maxtries:
                    self._report = [p1, trial, lamb]
                    return p1, jtj.diagonal(), message
                try:
                    jtjdiag = jtj.diagonal()
                    jtj[np.diag_indices_from(jtj)] += lamb
                    self.d0 = - solve(jtj, jtr)

                    if geo:
                        self._rpp = 2./h*((self.res(p1+h*self.d0) - r0)/h
                                          - j.dot(self.d0))
                        self.d1 = - solve(jtj, j.T.dot(self._rpp))/2
                        truncerr = 2*(np.sqrt(self.d1.dot(self.d1))/
                                      np.sqrt(self.d0.dot(self.d0)))

                    jtj[np.diag_indices_from(jtj)] = jtjdiag

                except LinAlgError as er:
                    print("\tSingular matrix, lamb = {}".format(lamb))
                    jtj[np.diag_indices_from(jtj)] -= lamb
                    lamb *= reject
                    tries += 1
                    continue

                lamb = np.nan_to_num(lamb)

                trial = p1 + self.d0
                if geo:
                    trial += self.d1

                r1 = self.res(trial, *self.args)
                nlnprob1 = r1.dot(r1)/2.
                 
                if geo and truncerr > 2:
                    pass
                elif nlnprob1 < nlnprob0: #success
                    success = True
                else: #failed step
                    pass

                if success:
                    lamb /= accept
                    p1 = trial
                else:
                    lamb *= reject

                tries += 1
                if iprint > 1:
                    print("\tsuccess = {}, lnprob = {}, \
                          lamb = {}".format(success, nlnprob1, lamb))

            #Convergence condition
            if (nlnprob0 - nlnprob1)/nlnprob0 < delta:
                if iprint:
                    print("log-prob relative changed by less than delta = {},\
                          optimum found".format(delta))
                message += 1
                break
        message += 1
        if iprint:
            print(self.msgdict[message])
        return p1, jtj.diagonal(), message

def lm_v2(x, y, pars, funct, ilambda = 10, niter = 20, njacobian = True,
    weights = 1.0, w_cut = 1e-10, chi2_stop = 1e-3, cvm = False,istep = 10.,
    fix = None,limits = None, autolambda = False, accel = False, **kwargs):
    """
    Levemberg Marquardt algorithm
    D. Orozco Suarez (summer 2021)
    All inputs numpy arrays please, be serious or use C.
    x - independent values
    y - dependent values
    pars - parameters (initial estimate)
    funct - funtion to fit
        called as y, jac = funct(x,pars) if njacobian = False
        called as y, _   = funct(x,pars) if njacobian = True
    njacobian = True ; estimate jacobian numerically
    weights = array of same size of x if given. == 1 if not given.
    w_cut = SVD inverse cut (whatever)
    ilambda = 10. ; initial lambda parameter. Could be 0.1 but who cares.
    niter = 20. ; max default iterations
    chi2_stop = 1e-3 ; stopping criteria
        when np.abs((ochi2 - chi2)/chi2)*100 < chi2_stop

    Example: lm_test() is just a full example.
    """
    if autolambda:
        print('Auto lambda')
    def __orderOfMagnitude(number):
        return 10**np.floor(np.log10(number))

    def check_limits(input,limits):
        if limits.any() != None:
            if len(limits.flatten()) > 3:
                for lim in np.arange(len(limits[:,0])):
                    pos = int(lim)
                    idx = int(limits[pos,0])
                    if input[idx] > limits[pos,2]:
                        input[idx] = limits[pos,2]
                    if input[idx] < limits[pos,1]:
                        input[idx] = limits[pos,1]
            else:
                idx = int(limits[0])
                if input[idx] > limits[2]:
                    input[idx] = limits[2]
                if input[idx] < limits[1]:
                    input[idx] = limits[1]

        return input

    if limits:
        check_limits(pars,limits)

    # calculate jacobian
    if njacobian:
        #numerically calculate the jacobian
        yfit, jac = numerical_der(x,pars,funct, **kwargs)
    else:
        yfit, jac = funct(x,pars, **kwargs)

    #Check length of x (can be anything provided it is flatten()
    x_length = len(yfit)
    pars_length = len(pars)

    if isinstance(weights, float) :
        print("weights is float")
        w = np.ones(x_length)
    else:
        w = np.copy(weights)

    print(fix)
    if fix is None:
        print("Fix is None")
        fix = np.ones(pars_length)
    else:
        print("Fix is NOT None")
        if len(fix) != pars_length:
            print('Fix ne x_length')
            return

    free = x_length - pars_length
    if free <= 1:
        print('not enough points')
        return

    #Set derivaties to zero when not taken into account
    jac = jac * fix[np.newaxis,:]
    #determine jacobian of merit function
    ochi = (y - yfit) * w
    J = np.matmul(ochi, jac)
    H = np.matmul(np.transpose(jac), jac*w[:,np.newaxis])
    ochi2 = np.sum(ochi**2)/free
    loop = 0
    if autolambda:
        ilambda = __orderOfMagnitude(autolambda*np.sqrt(np.linalg.norm(J)))

    while loop < niter:

        if cvm:
            covar = np.sqrt(np.outer(H.diagonal(),H.diagonal()))
            H /= covar
            H = np.nan_to_num(H)
            np.fill_diagonal(H, (1+ilambda))
            Hi = svd_solve(H,w_cut =w_cut) / covar
            Hi = np.nan_to_num(Hi)
            delta = np.matmul(Hi,J)
        else:
            np.fill_diagonal(H, H.diagonal()*(1+ilambda))
            delta = svd_solve(H, b = J,w_cut =w_cut)

        if accel:
            hpar = 1e-5
            # yfit_prev , _ = funct(x,pars, **kwargs) #f(x)
            yfit_accel, _ = funct(x,pars + hpar*delta*fix, **kwargs) #f(x+hd)
            chi_accel = (y - yfit_accel) * w             
            accel_val = 2./hpar * ( (chi_accel - ochi)/hpar - J.T.dot(delta) )
            accel_val = svd_solve(jac.T.dot(jac), b = jac.T.dot(accel_val),w_cut =w_cut)
            print(accel_val)

        new_pars = pars + delta * fix
        
        if limits:
            check_limits(new_pars,limits)

        yfit, _ = funct(x,new_pars, **kwargs)
        chi = (y - yfit) * w
        chi2 = np.sum(chi**2)/free

        if chi2 - ochi2 < 0:
            print('{:<6s}{:>3.0f}{:<8s}{:>12.4e}{:<6s}{:>1.2e}{:<6s}'.format('Iter: ',loop,' Lambda: ',ilambda,' chi2: ',ochi2,' better'))
            ilambda /= istep
            pars = np.copy(new_pars)

            # calculate jacobian
            if njacobian:
                #numerically calculate the jacobian
                yfit, jac = numerical_der(x,pars,funct, **kwargs)
            else:
                yfit, jac = funct(x,pars, **kwargs)
            jac = jac * fix[np.newaxis,:]


            #determine jacobian of merit function
            ochi = (y - yfit) * w
            J = np.matmul(ochi, jac)
            H = np.matmul(np.transpose(jac), jac*w[:,np.newaxis])
            if np.abs((ochi2 - chi2)/chi2)*100 < chi2_stop:
                print('STOP because (ochi2 - chi2)/chi2)*100 < chi2_stop')
                break
            ochi2 = np.sum(ochi**2)/free

        else:
            print('{:<6s}{:>3.0f}{:<8s}{:>12.4e}{:<6s}{:>1.2e}{:<6s}'.format('Iter: ',loop,' Lambda: ',ilambda,' chi2: ',ochi2,' worse'))
            ilambda *= istep

        if (ilambda < 1e-12) or (ilambda > 1e12):
            print('STOP because ilambda reached a limit')
            break
        loop += 1
    if loop == niter:
        print('STOP because max niter')

    chi2 = np.sum(chi**2)/free
    Hi = svd_solve(H)
    sigma = np.sqrt(Hi.diagonal())

    return pars, yfit, sigma, chi2

def lm(x, y, pars, funct, ilambda = 10, niter = 20, njacobian = True,
    weights = 1.0, w_cut = 1e-10, chi2_stop = 1e-3, cvm = False,istep = 10.,
    fix = None,limits = None, autolambda = False, accel = False, **kwargs):
    """
    Levemberg Marquardt algorithm
    D. Orozco Suarez (summer 2021)
    All inputs numpy arrays please, be serious or use C.
    x - independent values
    y - dependent values
    pars - parameters (initial estimate)
    funct - funtion to fit
        called as y, jac = funct(x,pars) if njacobian = False
        called as y, _   = funct(x,pars) if njacobian = True
    njacobian = True ; estimate jacobian numerically
    weights = array of same size of x if given. == 1 if not given.
    w_cut = SVD inverse cut (whatever)
    ilambda = 10. ; initial lambda parameter. Could be 0.1 but who cares.
    niter = 20. ; max default iterations
    chi2_stop = 1e-3 ; stopping criteria
        when np.abs((ochi2 - chi2)/chi2)*100 < chi2_stop

    Example: lm_test() is just a full example.
    """
    if autolambda:
        print('Auto lambda')
    def __orderOfMagnitude(number):
        return 10**np.floor(np.log10(number))

    def check_limits(input,limits):
        if limits.any() != None:
            if len(limits.flatten()) > 3:
                for lim in np.arange(len(limits[:,0])):
                    pos = int(lim)
                    idx = int(limits[pos,0])
                    if input[idx] > limits[pos,2]:
                        input[idx] = limits[pos,2]
                    if input[idx] < limits[pos,1]:
                        input[idx] = limits[pos,1]
            else:
                idx = int(limits[0])
                if input[idx] > limits[2]:
                    input[idx] = limits[2]
                if input[idx] < limits[1]:
                    input[idx] = limits[1]

        return input

    if limits:
        check_limits(pars,limits)

    # calculate jacobian
    if njacobian:
        #numerically calculate the jacobian
        yfit, jac = numerical_der(x,pars,funct, **kwargs)
    else:
        yfit, jac = funct(x,pars, **kwargs)

    #Check length of x (can be anything provided it is flatten()
    x_length = len(yfit)
    pars_length = len(pars)

    if isinstance(weights, float) :
        print("weights is float")
        w = np.ones(x_length)
    else:
        w = np.copy(weights)

    print(fix)
    if fix is None:
        print("Fix is None")
        fix = np.ones(pars_length)
    else:
        print("Fix is NOT None")
        if len(fix) != pars_length:
            print('Fix ne x_length')
            return

    free = x_length - pars_length
    if free <= 1:
        print('not enough points')
        return

    #Set derivaties to zero when not taken into account
    jac = jac * fix[np.newaxis,:]
    #determine jacobian of merit function
    ochi = (y - yfit) * w
    J = np.matmul(ochi, jac)
    H = np.matmul(np.transpose(jac), jac*w[:,np.newaxis])
    ochi2 = np.sum(ochi**2)/free
    loop = 0
    if autolambda:
        ilambda = __orderOfMagnitude(autolambda*np.sqrt(np.linalg.norm(J)))

    while loop < niter:

        if cvm:
            covar = np.sqrt(np.outer(H.diagonal(),H.diagonal()))
            H /= covar
            H = np.nan_to_num(H)
            np.fill_diagonal(H, (1+ilambda))
            Hi = svd_solve(H,w_cut =w_cut) / covar
            Hi = np.nan_to_num(Hi)
            delta = np.matmul(Hi,J)
        else:
            np.fill_diagonal(H, H.diagonal()*(1+ilambda))
            delta = svd_solve(H, b = J,w_cut =w_cut)

        new_pars = pars + delta * fix
        
        if limits:
            check_limits(new_pars,limits)

        yfit, _ = funct(x,new_pars, **kwargs)
        chi = (y - yfit) * w
        chi2 = np.sum(chi**2)/free

        if chi2 - ochi2 < 0:
            print('{:<6s}{:>3.0f}{:<8s}{:>12.4e}{:<6s}{:>1.2e}{:<6s}'.format('Iter: ',loop,' Lambda: ',ilambda,' chi2: ',ochi2,' better'))
            ilambda /= istep
            pars = np.copy(new_pars)

            # calculate jacobian
            if njacobian:
                #numerically calculate the jacobian
                yfit, jac = numerical_der(x,pars,funct, **kwargs)
            else:
                yfit, jac = funct(x,pars, **kwargs)
            jac = jac * fix[np.newaxis,:]


            #determine jacobian of merit function
            ochi = (y - yfit) * w
            J = np.matmul(ochi, jac)
            H = np.matmul(np.transpose(jac), jac*w[:,np.newaxis])
            if np.abs((ochi2 - chi2)/chi2)*100 < chi2_stop:
                print('STOP because (ochi2 - chi2)/chi2)*100 < chi2_stop')
                break
            ochi2 = np.sum(ochi**2)/free

        else:
            print('{:<6s}{:>3.0f}{:<8s}{:>12.4e}{:<6s}{:>1.2e}{:<6s}'.format('Iter: ',loop,' Lambda: ',ilambda,' chi2: ',ochi2,' worse'))
            ilambda *= istep

        if (ilambda < 1e-12) or (ilambda > 1e12):
            print('STOP because ilambda reached a limit')
            break
        loop += 1
    if loop == niter:
        print('STOP because max niter')

    chi2 = np.sum(chi**2)/free
    Hi = svd_solve(H)
    sigma = np.sqrt(Hi.diagonal())

    return pars, yfit, sigma, chi2

def numerical_der(x,pars,funct,**kwargs):

    #    for key, value in kwargs.items():
    try:
        h = kwargs['h']
    except Exception:
        h = 1
    y, _ = funct(x,pars, **kwargs)
    perturbation = np.copy(pars)
    y_length = len(y)
    pars_length = len(pars)
    jac = np.zeros((y_length,pars_length))

    for i in range(pars_length):
    #     if abs(pars[i]) > 1e-9:
    #         perturbation[i] = pars[i] * (1. + h)
    #         y_d, _ = funct(x,perturbation)
    #         perturbation[i] = pars[i] / (1. + h) * (1. - h)
    #         y_i, _ = funct(x,perturbation)
    #         perturbation[i] = pars[i] / (1. - h)
    #     else:
        perturbation[i] = pars[i] + h
        y_d, _ = funct(x,perturbation, **kwargs)
        perturbation[i] = pars[i] - 2*h
        y_i, _ = funct(x,perturbation, **kwargs)
        perturbation[i] = pars[i] + h
        jac[:,i] = (y_d - y_i)/(2.*h)

    return y,jac

def lm_test(cvm=True,istep=10.,niter=100,njacobian=True,autolambda=False,limits = None):
    def test_func(x,a, **kwargs):
        x_length = len(x)
        pars_length = len(pars)
        jac = np.zeros((x_length,pars_length))
        y = a[0] + a[1]*x + a[2]*x**2
        jac[:,0] = 1
        jac[:,1] = x
        jac[:,2] = x**2
        return y,jac
    x = np.arange(50)
    pars = [.75,1.70,2.45]
    y,jac = test_func(x,pars)
    y = y + y*np.random.uniform(-1,1,50)*1e-2
    pars = [.70,1.10,0.45]
    #new_pars,yfit = lm(x, y, pars, test_func)
    pars,yfit , _  ,_ = lm(x, y, pars,test_func,njacobian=njacobian,niter=niter,cvm=cvm,istep=istep,autolambda=autolambda,limits=limits)
    plt.plot(x,y,'o')
    plt.plot(x,yfit)
    print([.75,1.70,2.45],pars)
    plt.show()
    plt.plot(x,y-yfit,'.-')
    plt.show()

def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            t = (te - ts) * 1000
            if t > 1e3:
                print('%r  %2.2f seconds' %
                      (method.__name__, t/1000.))
            elif t > 1e6:
                print('%r  %2.2f mim' %
                      (method.__name__, t/1000./60.))
            else:
                print('%r  %2.2f ms' %
                      (method.__name__, t))

        return result
    return timed

def svd_solve(A, b=None, w_cut=1e-10):
    """
    This function solves the system of equations Ax=b by calculating the
    inverse of A using the SVD method: x=A^(-1)*b; A^(-1)=V*S^(-1)*U'
    Inputs:
        A: 2D array of dimensions nxm (n>=m)
        b: 1D array of dimensions n
        w_cut: cut-off frequency for singular values (fraction of the maximum).
        Diagonal elements S^(-1) are zero for the positions of S where its
        value is less than w_cut.
    """
    # Ab = np.abs(A)
    # A = np.where(Ab < np.min(Ab)*10.,0,A)
    # #Ac = np.zeros_like(b)
    # #Ac[:,0] = A[:,0]
    # #b = np.where(Ac == 0.0, 0 ,b)
    # bb = np.abs(b)
    # b = np.where(bb < np.min(bb)*10.,0,b)

    U, S, Vt = np.linalg.svd(A)
    sigma = w_cut*np.max(S)
    Sinv = np.where(S < sigma, 0, (1/S))
    # print(S)
    # print(sigma)
    zeros_Sinv = np.argwhere(Sinv == 0)
    Ainv = np.dot(np.transpose(Vt)*Sinv, np.transpose(U))
    return Ainv if b is None else np.dot(Ainv, b)

def cal_unit(pol_angle, retardance, ret_angle, angle_rot=None):
    """
    #TODO implement derivatives (it is in lib_pmp_v1.4)
    #theta = off set of the retarder
    #delta = retardance
    #alpha = off set of the polarizer(angle)
    """
    theta_r = ret_angle * np.pi/180.
    delta_r = retardance * np.pi/180.
    alpha_r = pol_angle * np.pi/180.
    c2 = np.cos(2*theta_r)
    s2 = np.sin(2*theta_r)

    return (
        0.5
        * np.matrix(
            [
                1.0,
                (c2**2.0 + s2**2 * np.cos(delta_r)) * np.cos(2 * alpha_r)
                + c2 * s2 * (1 - np.cos(delta_r)) * np.sin(2 * alpha_r),
                (s2**2.0 + c2**2 * np.cos(delta_r)) * np.sin(2 * alpha_r)
                + c2 * s2 * (1 - np.cos(delta_r)) * np.cos(2 * alpha_r),
                s2 * np.sin(delta_r) * np.cos(2 * alpha_r)
                - c2 * np.sin(delta_r) * np.sin(2 * alpha_r),
            ]
        ).T
    )

def test_cal_unit():
    theta = np.arange(0, 360, 1)
    out = np.zeros((4,360))
    for i in theta:
        out[:,i] = cal_unit(0.,75.,i,angle_rot=0).flatten()
    plt.plot(theta, out[0, :], label = 'I1')
    plt.plot(theta, out[1, :], label = 'I2')
    plt.plot(theta, out[2, :], label = 'I3')
    plt.plot(theta, out[3, :], label = 'I4')
    plt.title('rot = 0 deg')
    plt.legend()
    plt.show()

    for i in theta:
        out[:, i] = cal_unit(0., 75., i, angle_rot=60).flatten()
    plt.plot(theta, out[0, :], label='I1')
    plt.plot(theta, out[1, :], label='I2')
    plt.plot(theta, out[2, :], label='I3')
    plt.plot(theta, out[3, :], label='I4')
    plt.title('rot = 60 deg')
    plt.legend()
    plt.show()

def rotation_matrix(angle_rot):
    c, s = np.cos(2*angle_rot*np.pi/180), np.sin(2*angle_rot*np.pi/180)
    return np.matrix([[1, 0, 0, 0], [0, c, s, 0], [0, -s, c, 0], [0, 0, 0, 1]])

def test_rotation_matrix():
    return

def polarizer(alpha):
    #OJO ALPHA IS = 0 WHEN THE POLARIZER IS HORIZONTAL!!!! Vertical = > THETA = 90
    #this is a right handed system
    # y
    # |
    # |
    # |
    # --------- -> X

    cd = np.cos( np.deg2rad(2*alpha))
    sd = np.sin( np.deg2rad(2*alpha))
    c2 = cd**2
    s2 = sd**2

    return 0.5 * np.matrix(
        [
            [1, cd, sd, 0],  # --> row1 pl[0,:]
            [cd, c2, cd * sd, 0],  # --> row2 pl[1,:]
            [sd, cd * sd, s2, 0],  # --> row3 pl[2,:]
            [0, 0, 0, 0],
        ]
    )

def retarder(alpha,retardance):
    '''
    a retarder
    '''
    #OJO ALPHA IS = 0 fast axis horizontal
    #this is a right handed system
    # y
    # |
    # |
    # |
    # --------- -> X

    cd = np.cos( np.deg2rad(2*alpha) )
    sd = np.sin( np.deg2rad(2*alpha) )
    c2 = cd**2
    s2 = sd**2
    cr = np.cos( np.deg2rad(retardance))
    sr = np.sin( np.deg2rad(retardance))

    return np.matrix(
        [
            [1, 0, 0, 0],  # --> row1 rt[0,:]
            [0, c2 + s2 * cr, cd * sd * (1 - cr), -sd * sr,],  # --> row2 rt[1,:]
            [0, cd * sd * (1 - cr), s2 + c2 * cr, cd * sr],  # --> row3 rt[2,:]
            [0, sd * sr, -cd * sr, cr],
        ]
    )

def pol_lin(alpha):
    #OJO ALPHA IS = 0 WHEN THE POLARIZER IS HORIZONTAL!!!! Vertical = > THETA = 90
    #this is a right handed system
    # y
    # |
    # |
    # |
    # --------- -> X

    cd = np.cos( np.deg2rad(2*alpha))
    sd = np.sin( np.deg2rad(2*alpha))
    c2 = cd**2
    s2 = sd**2

    return 0.5 * np.matrix(
        [
            [1, cd, sd, 0],
            [cd, c2, cd * sd, 0],
            [sd, cd * sd, s2, 0],
            [0, 0, 0, 0],
        ]
    )

def rotating_wp(alpha,retardance):
    '''
    a retarder
    '''
    #OJO ALPHA IS = 0 fast axis horizontal
    #this is a right handed system
    # y
    # |
    # |
    # |
    # --------- -> X

    cd = np.cos( np.deg2rad(2*alpha))
    sd = np.sin( np.deg2rad(2*alpha))
    c2 = cd**2
    s2 = sd**2
    cr = np.cos( np.deg2rad(retardance))
    sr = np.sin( np.deg2rad(retardance))

    return np.matrix(
        [
            [1, 0, 0, 0],
            [0, c2 + s2 * cr, cd * sd * (1 - cr), -sd * sr],
            [0, cd * sd * (1 - cr), s2 + c2 * cr, cd * sr],
            [0, sd * sr, -cd * sr, cr],
        ]
    )

def test_rotating_wp(retardance=150.):

    theta = np.arange(0, 360, 1)
    polange = np.arange(0, 90, 0.25)
    out = np.zeros((360, 4))
    for i in range(len(theta)):
        out[i, :] = np.matmul(pol_lin(polange[0]), rotating_wp(
            theta[i], retardance))[0, :]
    plt.plot(theta, out[:, 0], label='I1')
    plt.plot(theta, out[:, 1], label='I2')
    plt.plot(theta, out[:, 2], label='I3')
    plt.plot(theta, out[:, 3], label='I4')
    #TODO ESTABA POR AQUI
    I2 = 0.5*(np.cos(2*theta*np.pi/180.)**2+np.sin(2*theta*np.pi /
                                               180.)**2*np.sin(retardance*np.pi/180))
    I3 = 0.5*np.cos(2*theta*np.pi/180.)*np.sin(2*theta*np.pi /
                                               180.)*(1-np.cos(retardance*np.pi/180))
    I4 = -0.5*np.sin(2*theta*np.pi/180.)*np.sin(retardance*np.pi/180)
    plt.plot(theta, I2, label='I2b', linestyle='--')
    plt.plot(theta, I3, label='I3b', linestyle='--')
    plt.plot(theta, I4, label='I4b', linestyle='--')
    plt.title(f'rot = {str(retardance)} deg')
    plt.legend()
    plt.show()

def lcvr(theta,delta):

    '''
    lcvr is just a retarder
    '''
    #quarter-wave plate, fast axis horizontal  THETA = 0, DELTA = 90
    #quarter-wave plate, fast axis vertical  THETA = 90, DELTA = 90
    #theta is fast axis, delta retardance
    #angles are measured positive counterclockwise
    #z axis positive direction is along the direction of light propagation
    #this is a right handed system
    # y
    # |
    # |
    # |
    # --------- -> X
    # angle in degree (changed internally into rad)
    # matrix [row,columns]

    deltar = np.deg2rad(delta)
    thetar = np.deg2rad(theta)
    c2 = np.cos(2*thetar)
    s2 = np.sin(2*thetar)

    return np.matrix([[1 , 0                         , 0                          , 0                  ],
                      [0 , c2**2+s2**2*np.cos(deltar), c2*s2*(1-np.cos(deltar))   , -s2*np.sin(deltar) ],
                      [0 , c2*s2*(1-np.cos(deltar))  , s2**2+c2**2*np.cos(deltar) , c2*np.sin(deltar)  ],
                      [0 , s2*np.sin(deltar)         , -c2*np.sin(deltar)         , np.cos(deltar)    ]])

def instrument_model(pardata):

    # get variables from global if not given
    delta1 = pardata['delta1']
    delta2 = pardata['delta2']
    theta1 = pardata['theta1']
    theta2 = pardata['theta2']
    pol_angle = pardata['pol_angle']
    rot_inst = pardata['rot_inst']

    mod_matrix = np.zeros((4,4),dtype=np.float64)
    for i in range(4):
        LC1 = retarder(theta1,delta1[i])
        LC2 = retarder(theta2,delta2[i])
        PL = polarizer(pol_angle)
        MR = rotation_matrix(rot_inst)
        dummy = np.matmul(np.matmul(PL, np.matmul(LC2, LC1)),MR)
        mod_matrix[i, :] = dummy[0,:]
    return mod_matrix

def lm_pol_model(theta,input_parameters,plot=False,**kwargs):

    '''
    #TODO implement derivatives (it is in lib_pmp_v1.4)
    #theta = off set of the retarder
    #delta = retardance
    #alpha = off set of the polarizer(angle)
    '''
    try:
        modulation = kwargs['modulation']
    except Exception:
        modulation = None

    pars = {
        'delta1': input_parameters[:4],
        'delta2': input_parameters[4:8],
        'theta1': input_parameters[8],
        'theta2': input_parameters[9],
        'pol_angle': input_parameters[10],
        'rot_inst': input_parameters[11],
    }
    transmission = input_parameters[12:16]

    alpha = input_parameters[16]
    delta = input_parameters[17]
    angle_rot = input_parameters[18]

    out = np.zeros((len(theta),4))
    pm = 2.*instrument_model(pars)

    for i in range(len(theta)):
        dummy = np.matmul(pm, cal_unit(
            alpha, delta, theta[i], angle_rot=angle_rot))
        out[i, :] = dummy.flatten()
    if plot:
        _extracted_from_pol_cal_model_(theta, out)
    jac = 0

    for i in range(4):
        out[:,i] *= transmission[i]

    if modulation != None:
        return out[:,modulation].flatten(order='F'), jac
    else:
        return out.flatten(order='F'), jac

def pol_cal_model(pars,alpha=0, delta=75, theta=0, angle_rot=None,plot=None):

    if not isinstance(theta, np.ndarray):
        theta=np.arange(0, 360, 1)
    '''
    #TODO implement derivatives (it is in lib_pmp_v1.4)
    #theta = off set of the retarder
    #delta = retardance
    #alpha = off set of the polarizer(angle)
    '''

    #    alpha = 0
    #    delta = 75
    #    theta = np.arange(0, 360, 1)
    #    polange = np.arange(0, 90, 0.25)
    out = np.zeros((len(theta),4))
    pm = 2.*instrument_model(pars)
    # ipm = svd_solve(pm)
    for i in range(len(theta)):
        dummy = np.matmul(pm, cal_unit(
            alpha, delta, theta[i], angle_rot=angle_rot))
        out[i, :] = dummy.flatten()
    if plot is not None:
        _extracted_from_pol_cal_model_(theta, out)
    return pm, out


# TODO Rename this here and in `pol_cal_model` and in `lm_pol_model`
def _extracted_from_pol_cal_model_(theta, out):
    plt.plot(theta, out[:, 0], label='I1')
    plt.plot(theta, out[:, 1], label='I2')
    plt.plot(theta, out[:, 2], label='I3')
    plt.plot(theta, out[:, 3], label='I4')
    plt.legend()
    plt.show()

@timeit
def random_model(sample=10000, theta = 0):
    if isinstance(theta,int):
        theta = np.arange(0, 360, 1)
        size_sample = 360
    else:
        size_sample = theta.shape
    arr = np.random.rand(12, sample)*360
    modeli = np.zeros((sample, size_sample, 4))
    modelo = np.zeros(( sample, 4, 4))

    for i in range(sample):
        pardata['delta1'] = arr[0:4, i]
        pardata['delta2'] = arr[4:8, i]
        pardata['rot_inst'] = arr[8, i]
        pardata['theta1'] = 0#arr[9, i]/360*1-0.5
        pardata['theta2'] = 45#arr[10, i]/360*1+44.5
        pardata['pol_angle'] = 90#arr[11, i]
        a, b = pol_cal_model(pardata, theta=theta)

        modeli[i,:,:] = b
        modelo[i,:,:] = a

    #np.random.shuffle(arr)
    return modeli, modelo

def mirror(theta,t):
    '''
    Capitani et al, 1989, Sol. Phys., 120, 173
    '''
    a = theta * np.pi/180.

    return np.array([[ (t**2. + 1) , (t**2. - 1) ,     0      ,     0      ],
        [  (t**2. - 1) , (t**2. + 1) ,     0      ,     0      ],
        [       0      ,      0      , 2*t*np.cos(a) , 2*t*np.sin(a) ],
        [       0      ,      0      , -2*t*np.sin(a), 2*t*np.cos(a) ]]) / 2.

def smirror(Phase_shift,R, Angle_incidence):
    '''
    '''
    #reflectance = r
    #Phase_shift
    alpha = Phase_shift * np.pi/180. #small angle app
    #angle of incidence
    beta = Angle_incidence * np.pi/180.

    a = (R+1.)/2.
    b = (R-1.)/2.

    return [[ a          , b*np.cos(2.*beta)                                                 ,    b*np.sin(2.*beta)                                              ,     0                                     ],
    [  b*np.cos(2.*beta) , a*np.cos(2.*beta)**2.-np.sqrt(R)*np.cos(alpha)*np.sin(2.*beta)**2 , (a+np.sqrt(R)*np.cos(alpha))*np.sin(4.*beta)/2.                   ,  np.sqrt(R)*np.sin(alpha)*np.sin(2.*beta) ],
    [  b*np.cos(2.*beta) , (a+np.sqrt(R)*np.cos(alpha))*np.sin(4.*beta)/2.                   , a*np.sin(2.*beta)^2.-np.sqrt(R)*np.cos(alpha)*np.cos(2.*beta)**2. , -np.sqrt(R)*np.sin(alpha)*np.cos(2.*beta) ],
    [       0            , -np.sqrt(R)*np.sin(alpha)*np.sin(2.*beta)                         , np.sqrt(R)*np.sin(alpha)*np.cos(2.*beta)                          , -np.sqrt(R)*np.cos(alpha)                 ]]


def fit_calib(x,y,plot=None,w_cut=1e-10):
    '''
    input:  calib states (angle,4), obs (angle)
    '''
    xc = np.matmul(x.T,x)
    xi = svd_solve(xc, w_cut=w_cut)
    cf = np.matmul(xi,np.matmul(x.T,y))
    yfit = np.matmul(x,cf)
    cov = np.sqrt(np.sum((yfit - y)**2)/(len(y) - len(xc) + 1.0 ) * xc.diagonal() )
    if plot:
        plt.plot(y,'.')
        plt.plot(yfit,'-')
    return cf,yfit,cov

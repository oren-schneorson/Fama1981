
# Python script file fminunc.py:

import sys
#import numba as nb
import numpy as np
from numpy.linalg import eig, det, inv
from numpy.linalg import cholesky as chol
from scipy.sparse import spdiags
from scipy import optimize
#from python_functions import *


#algo = sys.argv[1]

x = np.loadtxt('X.csv', delimiter=',')
y = np.loadtxt('Y.csv', delimiter=',')

beta0 = np.loadtxt('beta0.csv', delimiter=',')
N = len(y)
y = np.reshape(y, (N, 1))


# init conditions
theta_prev = -1
theta_ML = -.9
precision = 1e-5
precision_theta = 1e-2
precision_beta = 1e-2
spectrum_theta = [-1, -.2]
spectrum_beta = [-2, 2]

ub_theta = -0.1
lb_theta = -.999
ub_beta = 2
lb_beta = -2
bnds_theta = ((lb_theta, ub_theta), )
bnds_beta = ((lb_beta, ub_beta), )

beta_ML = np.array([-1.44, -0.23])



def get_Omega_movavg1(theta, N):
    # get Omegga with sigma2=1
    Omega = (1 + theta ** 2 - theta) * np.eye(N) + theta * np.ones((N, N))
    offsets = 3
    Omega1 = Omega - spdiags(theta*np.ones((2*(N-1)+1-offsets, N)), [i for i in range(-N+1, N) if abs(i) > (offsets-1)/2]).toarray()
    
    return Omega1


#@nb.njit(cache=False)
def get_Omega_movavg(Omega1, sigma2):
    
    eigen, _ = eig(Omega1)
    eigen = np.real(eigen)
    det_Omega = np.prod(np.multiply(eigen, sigma2))
    Omega = Omega1 * sigma2
    return Omega, det_Omega

#@nb.njit(cache=False)
def get_H_movavg_python(theta, N):
    '''
    GET_H_MOVAVG Get transformation of moving average model given theta.
       The Inverse of a Matrix Occurring in First-Order Moving-Average Models
       Relevant for any model with a residual of:
       epsilon_{t}+theta*epsilon_{t-1}, 
       where Cov(epsilon_{t}, epsilon_{t-1}) = 0
    
       For further info, see V. R. R. Uppuluri and J. A. Carpenter Sankhyā:
       The Indian Journal of Statistics, Series A (1961-2002)
       Mar., 1969, Vol. 31, No. 1 (Mar., 1969), pp. 79-82
       [Uppuluri-InverseMatrixOccurring-1969.pdf]
       Assumption: unit variance.
    '''
    

    if np.abs(theta) < 1e-15:
        # TODO: change this to nan
        V = np.empty((N,N))
        H = np.empty((N,N))
        return H, V

    #delta = np.array([1.0, 1 + theta ** 2] + [1.0]*(N-1)).astype(np.float64)

    #delta = [1, 1 + theta ** 2]
    #delta[1] = np.array([1 + theta ** 2]).astype(np.float64)
    #delta = np.ones((N+1, 1)).astype(np.float64)
    delta = np.ones((N+1, )).astype(np.float64)

    delta[1] = np.array(1 + theta ** 2).astype(np.float64)
    #delta[1] = 1.0 + theta ** 2

    for j in range(2, N+1):
        delta[j] = (1 + theta ** 2) * delta[j-1] - (theta ** 2) * delta[j-2]
        #delta.append( (1 + theta ** 2) * delta[-1] - (theta ** 2) * delta[-2] )

    v_1j = np.empty((N, ))
    v_jj = np.empty((N, ))
    #v_1j = list()
    #v_jj = list()

    for j in range(N):
        v_1j[j] = ((-1) ** j) * (theta ** j) * delta[N-1-j]/delta[N]
        v_jj[j] = delta[N-1-j] * delta[j]/delta[N]
        #v_1j.append( ((-1) ** j) * (theta ** j) * delta[N-1-j]/delta[N] )
        #v_jj.append( delta[N-1-j] * delta[j]/delta[N] )


    #v_1j = np.array(v_1j)
    #v_jj = np.array(v_jj)

    lambda_1j = np.divide(v_jj, v_1j)

    V = np.zeros((N, N))

    for j in range(1, N):
        V[j, j] = v_jj[j]
        for k in range(j+1,N):
            V[j, k] = lambda_1j[j] * v_1j[k]
            V[k, j] = V[j, k]

    if np.isnan(V).any() or np.isinf(V).any():
        #warn('Moving average parameter and sample size imply some elements could not be calculated due to computation error.')
        #warn('Setting matrix as symmetrical to the previous block up the diagonal of the matrix')

        V_bads = np.where(np.isinf(V) | np.isnan(V))
        #V_nans = np.where(np.isnan(V))
        #V_bads = V_infs | V_nans
        #V_bads = tuple(np.hstack((V_infs, V_nans)))

        V_computational_limit = V[min(V_bads[0]):, min(V_bads[1]):]
        V[min(V_bads[0]):, min(V_bads[1]):] = V[min(V_bads[0])-V_computational_limit.shape[0]:min(V_bads[0]), min(V_bads[1])-V_computational_limit.shape[1]:min(V_bads[1])]

    V[0, :] = np.reshape(v_1j, (N,))
    V[:, 0] = v_1j.transpose()
    H = chol(V)

    return H, V


#@nb.njit(cache=False)
def Omega_inv_hat__python(theta, N):
    H, _ = get_H_movavg_python(theta, N)
    return H.transpose() @ H


#@nb.njit(cache=False)
def Omega_inv_hat_python(theta, X, Y, beta, H, FGLS=False):
    N = X.shape[0]
    #Omega_inv_hat = Omega_inv_hat__python(theta, N)
    Omega_inv_hat = H.transpose() @ H
    if FGLS:
        beta = beta_FGLS_python(theta, X, Y, H)
    return Omega_inv_hat / sigma2_python(beta, X, Y)


#@nb.njit(cache=False)
def beta_FGLS_python(theta, X, Y, H):
    N = X.shape[0]
    K = X.shape[1]
    #Omega_inv_hat = Omega_inv_hat__python(theta, N)
    Omega_inv_hat = H.transpose() @ H
    return (inv(((X.transpose() @ Omega_inv_hat) @ X)) @ ((X.transpose() @ Omega_inv_hat) @ Y))

#@nb.njit(cache=False)
def sigma2_python(beta, X, Y):
    N = X.shape[0]
    err = Y - X@beta
    return np.sum(np.power(err, 2))/N


#@nb.njit(cache=False)
def L_python(theta, beta, X, Y, Omega1):

    beta = beta.reshape((len(beta), 1))
    N = X.shape[0]

    sigma2 = sigma2_python(beta, X, Y)
    H_hat, _ = get_H_movavg_python(theta, N)

    #Omega_inv_hat = Omega_inv_hat_python(theta, X, Y, beta, H)
    Omega_inv_hat = H_hat.transpose() @ H_hat
    Omega_inv_hat = Omega_inv_hat / sigma2

    err = Y - X@beta

    Omega_hat, _ = get_Omega_movavg(Omega1, sigma2)
    eigen, _ = eig(Omega_hat)
    eigen = np.real(eigen)

    output = -N/2*np.log(2*np.pi)-0.5*np.sum(np.log(eigen))-((err.transpose() @ Omega_inv_hat) @ err)/2/sigma2

    return output












Omega1 = get_Omega_movavg1(theta_ML, N)



# TNC: truncated Newton algorithm
# L-BFGS-B: L (local?) - Broyden-Fletcher-Goldfarb-Shanno - B (bounded?)
#algo = 'TNC'
#algo = 'L-BFGS-B'
while abs(theta_ML-theta_prev) > precision:
	
	theta_prev = theta_ML
	Omega1 = get_Omega_movavg1(theta_ML, N)

	fun_beta = lambda beta: -L_python(theta_ML, beta, x, y, Omega1)
	res_beta = optimize.minimize(fun_beta, beta0, method=algo, bounds=bnds_beta, options={'disp': True})
	beta_ML = res_beta.x
	

	fun_theta = lambda theta_: -L_python(theta_[0], beta_ML, x, y, Omega1)
	res_theta = optimize.minimize(fun_theta, theta_ML, method=algo, bounds=bnds_theta, options={'disp': True})
	theta_ML = res_theta.x[0]
	print('beta: ', beta_ML, 'theta:', theta_ML)
	print('*************')



matlab_output = [theta_ML, beta_ML[0], beta_ML[1]]
# beta:  [-1.44069486 -0.23299893] theta: -0.8525364
# beta:  [-1.44022318 -0.23319829] theta: -0.854433432631684

# result with including ILS/USD exchange rate didn't change the result of interest rate coefficient.

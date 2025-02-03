
import numpy as np
from scipy import optimize
from numpy.linalg import eig, det, inv
from numpy.linalg import cholesky as chol
from scipy.sparse import spdiags
from warnings import warn
import numba as nb



def get_Omega_movavg1(theta, N):
    # get Omegga with sigma2=1
    Omega = (1 + theta ** 2 - theta) * np.eye(N) + theta * np.ones((N, N))
    offsets = 3
    Omega1 = Omega - spdiags(theta*np.ones((2*(N-1)+1-offsets, N)), [i for i in range(-N+1, N) if abs(i) > (offsets-1)/2]).toarray()
    
    return Omega1



@nb.jit
def get_delta(theta, N):
    delta = [theta*0 + 1, 1 + theta ** 2]
    for j in range(2, N+1):
        delta.append( (1 + theta ** 2) * delta[-1] - (theta ** 2) * delta[-2] )
    return delta

@nb.jit
def get_vs(theta, N, delta):
    v_1j = list()
    v_jj = list()

    for j in range(N):
        v_1j.append( ((-1) ** j) * (theta ** j) * delta[N-1-j]/delta[N] )
        v_jj.append( delta[N-1-j] * delta[j]/delta[N] )
    return v_1j, v_jj



@nb.jit
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

    delta = get_delta(theta, N)

    v_1j, v_jj = get_vs(theta, N, delta)

    v_1j = np.array(v_1j).astype(np.float64)
    v_jj = np.array(v_jj).astype(np.float64)

    lambda_1j = np.divide(v_jj, v_1j)
    #lambda_1j = [v_jj[i] / v_1j[i] for i in range(len(v_jj))]

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


x = np.loadtxt('X.csv', delimiter=',')
y = np.loadtxt('Y.csv', delimiter=',')

N = x.shape[0]
K = x.shape[1]

beta_ML = np.array([-1.85450005, -0.23292365]).transpose()
Omega1 = get_Omega_movavg1(-.8, N)


@nb.jit
def fun_theta(theta):

    delta = np.ones((N+1, )).astype(np.float64)
    delta[1] = np.array(1 + theta[0] ** 2).astype(np.float64)
    for j in range(2, N+1):
        delta[j] = (1 + theta[0] ** 2) * delta[j-1] - (theta[0] ** 2) * delta[j-2]

    v_1j = np.empty((N, ))
    v_jj = np.empty((N, ))

    for j in range(N):
        v_1j[j] = ((-1) ** j) * (theta[0] ** j) * delta[N-1-j]/delta[N]
        v_jj[j] = delta[N-1-j] * delta[j]/delta[N]

    lambda_1j = np.divide(v_jj, v_1j)

    V = np.zeros((N, N))

    for j in range(1, N):
        V[j, j] = v_jj[j]
        for k in range(j+1,N):
            V[j, k] = lambda_1j[j] * v_1j[k]
            V[k, j] = V[j, k]

    V[0, :] = np.reshape(v_1j, (N,))
    V[:, 0] = v_1j.transpose()
    H_hat = chol(V)
    err = y - x@beta_ML
    sigma2 = np.sum(np.power(err, 2))/N

    Omega_inv_hat = H_hat.transpose() @ H_hat
    Omega_inv_hat = Omega_inv_hat / sigma2

    Omega_hat = Omega1 * sigma2

    eigen, _ = eig(Omega_hat)
    eigen = np.real(eigen)

    output = -N/2*np.log(2*np.pi)-0.5*np.sum(np.log(eigen))-((err.transpose() @ Omega_inv_hat) @ err)/2/sigma2
    print(theta, output)
    if np.isnan(output):
        raise ValueError

    return -output
    #return (theta+.8) ** 2


ub = -0.1
lb = -.999
bnds = ((lb, ub), )
res_theta = optimize.minimize(fun_theta, -.8, method='TNC', bounds=bnds, options={'disp': True})


#print(beta0_.shape)
#print(beta1_.shape)
print(res_theta)

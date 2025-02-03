
# Python script file fminunc.py:

import sys
from scipy import optimize
import numpy as np
from python_functions import *


#x = sys.argv[1]
#y = sys.argv[2]
#beta0 = sys.argv[3]
#beta0 = np.ndarray(beta0)

#x0 = np.reshape(x0, (len(x0), 1))
#x1 = np.reshape(x0, (len(x1), 1))

#x = np.concatenate([x0, x1], axis=1)

x = np.loadtxt('X.csv', delimiter=',')
y = np.loadtxt('Y.csv', delimiter=',')

if len(x.shape) == 1:
	x = x.reshape((x.shape[0], 1))

N = len(y)
y = np.reshape(y, (N, 1))


# init conditions
theta_prev = -1
theta_ML = -.9
precision = 1e-5


#beta_ML = np.array([-1.85450005, -0.23292365])
#beta_ML = np.array([-1, -0.23292365])
beta_ML = np.array([-1, ])
Omega1 = get_Omega_movavg1(theta_ML, N)

"""
@nb.jit
def fun_theta(theta_):
	return -L_python(theta_, beta_ML, x, y, Omega1)
"""


ub_theta = -0.1
lb_theta = -.999
ub_beta = 2
lb_beta = -2
bnds_theta = ((lb_theta, ub_theta), )
bnds_beta = ((lb_beta, ub_beta), )


# TNC: truncated Newton algorithm
# L-BFGS-B: L (local?) - Broyden-Fletcher-Goldfarb-Shanno - B (bounded?)
algo = 'TNC'
# algo = 'L-BFGS-B'
while abs(theta_ML-theta_prev) > precision:
	
	theta_prev = theta_ML
	Omega1 = get_Omega_movavg1(theta_ML, N)

	fun_beta = lambda beta: -L_python(theta_ML, beta, x, y, Omega1)
	print(fun_beta(beta_ML))
	res_beta = optimize.minimize(fun_beta, beta_ML, method=algo, bounds=bnds_beta, options={'disp': True})
	beta_ML = res_beta.x
	

	fun_theta = lambda theta_: -L_python(theta_, beta_ML, x, y, Omega1)
	res_theta = optimize.minimize(fun_theta, theta_ML, method=algo, bounds=bnds_theta, options={'disp': True})
	theta_ML = res_theta.x[0]
	print('beta: ', beta_ML, 'theta:', theta_ML)
	print('*************')


# beta:  [-1.44069486 -0.23299893] theta: -0.8525364
# beta:  [-1.44022318 -0.23319829] theta: -0.854433432631684

# result with including ILS/USD exchange rate didn't change the result of interest rate coefficient.

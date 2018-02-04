#!/usr/bin/python3

import sys
import linoptlib
import numpy as np
from scipy import optimize
#from sdaopt import sda

# number of input/output modes of the linear optical interfermoter
num_modes = 10

def print_fun(x, f, accepted):
	print("at minimum %.4f accepted %d" % (f, int(accepted)))

# optimize interferometer's unitary using specified method
if sys.argv[1] == 'BFGS':
	opt_val = optimize.minimize(linoptlib.cost_function, 6*np.random.rand(num_modes**2) - 3, method='BFGS', jac=False, tol=None, callback=None, options={'gtol':1e-5})
elif sys.argv[1] == 'Nelder-Mead':
	opt_val = optimize.minimize(linoptlib.cost_function, 6*np.random.rand(num_modes**2) - 3, method='Nelder-Mead', tol=1e-6, options={'maxfev':1e5})
elif sys.argv[1] == 'DiffEvol':
	bounds = np.tile(np.array([0,1]),(num_modes**2,1))
	opt_val = optimize.differential_evolution(linoptlib.cost_function, bounds, strategy='best1bin', maxiter=1000, popsize=10, tol=0.01, mutation=(0.5, 1), recombination=0.7, seed=None, callback=None, disp=False, polish=True, init='latinhypercube', atol=0)
elif sys.argv[1] == 'BasinHop':
	opt_val = optimize.basinhopping(linoptlib.cost_function, 6*np.random.rand(num_modes**2) - 3, niter=1000, T=0.02, stepsize=0.5, minimizer_kwargs={"method":"BFGS", "options":{'gtol':1e-5, 'norm':1}}, take_step=None, accept_test=None, callback=print_fun, interval=50, disp=False, niter_success=None)
elif sys.argv[1] == 'SDA':
    bounds = np.tile(np.array([0,1]),(num_modes**2,1))
    opt_val = sda(linoptlib.cost_function, 6*np.random.rand(num_modes**2) - 3, bounds, maxiter=1000, minimizer_kwargs={"method":"BFGS"}, initial_temp=10, visit=2.62, accept=-10.0, maxfun=1e7, seed=None, pure_sa=False)
else:
	opt_val = 'Wrong method'

# print result
#if isinstance(opt_val,str):
#	print(opt_val)
#else:
#	print(opt_val.fun)

#print(opt_val)
print('%s' %opt_val.fun + ' ' + ' '.join('%s' %x for x in opt_val.x))

# release memory allocated for stanisic_functor object
linoptlib.sf_destructor(linoptlib.sf)

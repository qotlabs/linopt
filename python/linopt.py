#!/usr/bin/python3

import ctypes
import numpy as np
from scipy import optimize
import math

linopt = ctypes.CDLL('../wrappers/linopt-wrap.so')

linopt.stanisic_functor_constructor.argtypes = (ctypes.c_int, ctypes.c_int, # Number of full system and ancilla modes
											ctypes.POINTER(ctypes.c_int), # Array of fock states, representing basis of the full system, Expected size = fb_size*Mf
											ctypes.c_int, # Number of fock states in the full basis
											ctypes.POINTER(ctypes.c_int), # Array of fock states, representing basis of the ancilla, Expected size = ab_size*Ma
											ctypes.c_int, # Number of fock states in the ancilla basis
											ctypes.POINTER(ctypes.c_int), # Input fock state. Expected size = Mf
											ctypes.POINTER(ctypes.c_int), # Two dimensional array of fock states. First index - target state index,
																			# second - fock state index in the given target state.
																			# Expected size = ts_size1*ts_size2*(Mf-Ma)
											ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), # Arrays of target states amlitudes (real and imaginary part)
																											# Expected size = ts_size1*ts_size2
											ctypes.c_int, ctypes.c_int) #Number of target states and nonzero amplitudes in target states respectively
linopt.stanisic_functor_constructor.restype = ctypes.c_void_p

sf_constructor = linopt.stanisic_functor_simple_constructor
sf_constructor.restype = ctypes.c_void_p

sf_destructor = linopt.stanisic_functor_destructor
sf_destructor.argtype = ctypes.c_void_p

costfunc = linopt.stanisic_functor_apply
costfunc.argtypes = (ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_int)
costfunc.restype = ctypes.c_double

full_photon_number = 4 # number of photons evolving in the linear-optical network
anc_photon_number = 2 # number of photons detected in the ancilla set of modes

Mf = 8 # full number of modes
Ma = 4 # ancilla number of modes

full_dim = math.factorial(Mf+full_photon_number-1)/(math.factorial(full_photon_number)*math.factorial(Mf-1)) # full fock state dimension 
anc_dim =math.factorial(Ma+anc_photon_number-1)/(math.factorial(anc_photon_number)*math.factorial(Ma-1)) # ancilla set dimension

ifs_array = np.array([1, 1, 1, 1, 0, 0, 0, 0]) # input fock state array
ifs_array_ptr = ifs_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)) # pointer to input fock state array

x_array = np.random.rand(64)
#x_array_ptr = x_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
x_size = x_array.size
x_array_ptr = (ctypes.c_double * x_size)(*x_array)

sf = sf_constructor()
res = costfunc(sf,x_array_ptr,x_size)

def cost_function(x):
	ar_size = x.size
	ar_ptr = (ctypes.c_double * ar_size)(*x)
	return -costfunc(sf,ar_ptr,ar_size)

#BFGS
opt_val = optimize.minimize(cost_function, np.random.rand(64), method='BFGS', jac=False, tol=1e-3, callback=None, options=None)
#Nelder-Mead
#opt_val = optimize.minimize(cost_function, np.random.rand(64), method='Nelder-Mead',tol=1e-4)
#differential evolution
a = np.array([0,1])
bounds = np.tile(a,(x_size,1))
#opt_val = optimize.differential_evolution(cost_function, bounds, strategy='best1bin', maxiter=1000, popsize=10, tol=0.01, mutation=(0.5, 1), recombination=0.7, seed=None, callback=None, disp=False, polish=True, init='latinhypercube', atol=0)

print(opt_val.fun)

sf_destructor(sf)

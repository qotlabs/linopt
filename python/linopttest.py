#!/usr/bin/python3

import sys
sys.path.insert(0, '../wrappers')
import random
import numpy as np
import linopttools as lot

from math import sqrt

#from minieigen import *
from pylinopt import *

from scipy.optimize import minimize

logical_basis, ancilla_basis = basis().generate_basis(2,4), basis().generate_basis(2,4)

full_basis = ancilla_basis*logical_basis

in_state = fock()
in_state[:] = [1,1,1,1,0,0,0,0]

# initialize fock states forming the bell basis
# these are used as keys to initialize the target state
fock_list    = [fock() for i in range(6)]
fock_list[0][:] = [1, 0, 1, 0]
fock_list[1][:] = [0, 1, 0, 1]
fock_list[2][:] = [1, 0, 0, 1]
fock_list[3][:] = [0, 1, 1, 0]
fock_list[4][:] = [1, 1, 0, 0]
fock_list[5][:] = [0, 0, 1, 1]

# initialize the list of dual-rail bell
# the target state
state_list = [state() for i in range(6)]
state_list[0][fock_list[0]] =  1/sqrt(2)
state_list[0][fock_list[1]] =  1/sqrt(2)
state_list[1][fock_list[0]] =  1/sqrt(2)
state_list[1][fock_list[1]] = -1/sqrt(2)
state_list[2][fock_list[2]] =  1/sqrt(2)
state_list[2][fock_list[3]] =  1/sqrt(2)
state_list[3][fock_list[2]] =  1/sqrt(2)
state_list[3][fock_list[3]] = -1/sqrt(2)
state_list[4][fock_list[4]] =  1/sqrt(2)
state_list[4][fock_list[5]] =  1/sqrt(2)
state_list[5][fock_list[4]] =  1/sqrt(2)
state_list[5][fock_list[5]] = -1/sqrt(2)

r = random.sample(range(100), 64)
r = [float(x)/100 for x in r]

c = circuit()
#u = unitary_matrix()

#a = lot.cf_inner_product(r, in_state, state_list, ancilla_basis)
a = minimize(lot.cf_inner_product, r, args = (in_state, state_list, ancilla_basis, full_basis, c), \
    method='BFGS', tol=None, options={'gtol':1e-5})
print(a)

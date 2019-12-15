#!/usr/bin/python3

import sys
sys.path.insert(0, '../pylib')
from linopt import *
import numpy as np
import tools as lot
from math import sqrt

from scipy.optimize import minimize

logical_basis = Basis(2, 4)
ancilla_basis = Basis(2, 4)

full_basis = ancilla_basis*logical_basis

in_state = Fock([1,1,1,1,0,0,0,0])

# initialize the list of dual-rail bell
# the target state
state_list = [
State({(1,0,1,0): 1/sqrt(2), (0,1,0,1):  1/sqrt(2)}),
State({(1,0,1,0): 1/sqrt(2), (0,1,0,1): -1/sqrt(2)}),
State({(1,0,0,1): 1/sqrt(2), (0,1,1,0):  1/sqrt(2)}),
State({(1,0,0,1): 1/sqrt(2), (0,1,1,0): -1/sqrt(2)}),
State({(1,1,0,0): 1/sqrt(2), (0,0,1,1):  1/sqrt(2)}),
State({(1,1,0,0): 1/sqrt(2), (0,0,1,1): -1/sqrt(2)}),
]

r = np.random.rand(64)

a = minimize(lot.cf_inner_product, r, args = (in_state, state_list, ancilla_basis, full_basis), \
    method='BFGS', tol=None, options={'gtol':1e-5})
print(a)

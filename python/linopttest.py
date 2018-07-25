import sys
sys.path.insert(0, '/home/dave/c++/linopt/wrappers')

from math import sqrt

from minieigen import *
from pylinopt import *

f = fock()
g = fock()
# logical_basis = basis()
# logical_basis.generate_basis(2,4,f)
ancilla_basis = basis()
ancilla_basis.generate_basis(2,4,g)
full_basis = basis()
full_basis.generate_basis(4,8,f)
print(len(full_basis))

in_state = fock()
in_state[:] = [1,1,1,1,0,0,0,0]

# initialize fock states forming the bell basis
fock_list    = [fock() for i in range(6)]
fock_list[0][:] = [1, 0, 1, 0]
fock_list[1][:] = [0, 1, 0, 1]
fock_list[2][:] = [1, 0, 0, 1]
fock_list[3][:] = [0, 1, 1, 0]
fock_list[4][:] = [1, 1, 0, 0]
fock_list[5][:] = [0, 0, 1, 1]

# initialize the list of dual-rail bell
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

v = VectorX.Random(64)
u = unitary_matrix()

c = chip()
c.input_state = in_state
c.output_basis = full_basis
c.unitary = u.exp_hermite(v)
print(len(c.output_basis))
postselected = state()

p, res = 0, 0

for anc in ancilla_basis:
    postselected = c.output_state().postselect(anc)
    p = postselected.norm()
    if p == 0:
        continue
    postselected /= p
    p = p*p
    for target_state in state_list:
        res += p * abs(postselected.dot(target_state))**10
    print(res)
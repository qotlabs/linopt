import random
from pylinopt import *
from math import sqrt

def cf_inner_product(random_vector, input_fock, target, ancilla, full_basis, c):
	"""
	The function computes the figure of merit of the output state optimization
	routine accounting both for fidelity of the output state and the target
	state and the probability of success of registering the desired state on the 
	output heralded by some ancilla subsystem state.

	Args:
		random_vector:  a 1xN^2 list of float numbers in [0,1] used to initialize
						a random starting point in the unitary matrix space

		input_fock:     input fock state, should be an instance of the class "fock"

		target:         the target state you want to achieve given an input fock state,
						a linear optical transformation and a list of ancilla fock
						states used to postselect the outcomes. Should be either an instance
		                of the class "state" or a list of instances of the class "state"

		ancilla:        an instance of the class "basis" containing a list of ancilla
						subsystem fock states which are used to postselect the output
						state in the logical subsystem. The number of photons and the
						number of modes in ancilla fock states should be consistent with
						the total number of modes and photons in the system and the number
						of photons and the modes restricted to the logical subsystem in
						the target argument

		full_basis:	    a list of fock states used to compute the output state before
						the postselection procedure. It may either contain the full fock state
						basis of the system or be composed only of the combination of
						the possible states in the logical and ancilla subsystems. Such a combination
						eliminates the necessity to compute the probabilities for the output
						states which won't be present in the system. For example, if the total
						number of photons in the system is 4, the ancilla subsystem is set to
						have 2 photons, and the logical subsystem is also set to have 2 photons,
						the states, where 3 or 4 photons populate either ancilla or logical subsystem
						are irrelevant to the problem being solved

	Returns:
		A float value computed accroding to the formula () accounting both for the fidelity
		of the linear optical transformation and the success probability.
	"""
	

	# check, if the user supplied non-empty arguments
	# if not random_vector:
	# 	raise ValueError('You supplied an empty vector to initialize the unitary matrix of the circuit')
	# if not input_fock:
	# 	raise ValueError('You supplied an empty fock input fock state')
	# if not target:
	# 	raise ValueError('You have supplied an empty target state')
	# if not ancilla:
	# 	raise ValueError('You have supplied an empty ancilla basis')

	# check the consistency between the vector size and the number of modes in input fock state
	# if len(input_fock.as_list) != sqrt(len(random_vector)):
	# 	raise ValueError('The supplied vector length is inconsistent with the number of input fock state modes')
	
	# check if the fock input state, output target state and ancilla basis are consistent ->
	# the total photon number in the input fock state equals total photon number of the target state +
	# + the total photon number in any of the ancilla basis states
	# for anc in ancilla:
	# 	for sta in target:
	# 		for key in sta.as_dict.keys():
	# 			if input_fock.total() != (anc.total() + key.total()):
	# 				raise ValueError('Input fock state total photon number is inconsistent with the target state and ancilla basis')
	# 			if len(input_fock.as_list) != (len(anc.as_list) + len(key.as_list)):
	# 				raise ValueError('Number of input modes is inconsistent with total output number of modes')   

	#v = VectorX(random_vector)
	c.circuit()
	c.input_state = input_fock
	c.output_basis = full_basis
	c.unitary = exp_hermite(random_vector)
	postselected = state()
	out_state = c.output_state()

	p, res = 0, 0

	out_postselect = out_state.postselect

	for anc in ancilla.as_list:
		postselected = out_postselect(anc)
		p = postselected.norm()
		if p == 0:
			continue
		postselected /= p
		p = p*p
		for target_state in target:
			res += p * abs(postselected.dot(target_state))**10

	return -res

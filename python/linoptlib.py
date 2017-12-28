import ctypes

# load the library containing computational tools for linear optics 
linopt = ctypes.CDLL('../wrappers/linopt-wrap.so')

# low-level stanisic_functor constructor
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

# stanisic_functor constructor with hard-coded initialization parameters
# Mf = 8, Ma = 4, full number of photons - 4, target states - 6 dual-rail Bell states
sf_constructor = linopt.stanisic_functor_simple_constructor
sf_constructor.restype = ctypes.c_void_p

# stanisic_functor destructor
sf_destructor = linopt.stanisic_functor_destructor
sf_destructor.argtype = ctypes.c_void_p

# cost function implementation based on the method described in Stanisic et al. Generating entanglement with linear optics, PRA 96, 043861 (2017)
# cost function return the probability of generation of any subset of 6 dual-rail Bell states set (two qubits - 4 modes and two photons) 
# given the configuration of the linear-optical interferometer (a nxn unitary matrix)
costfunc = linopt.stanisic_functor_apply
costfunc.argtypes = (ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_int)
costfunc.restype = ctypes.c_double

sf = sf_constructor()

# redefinition of the cost-function
def cost_function(x):
    ar_size = x.size
    ar_ptr = (ctypes.c_double * ar_size)(*x)
    return -costfunc(sf,ar_ptr,ar_size)
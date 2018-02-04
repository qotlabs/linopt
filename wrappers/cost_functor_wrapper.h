#ifndef COST_FUNCTOR_WRAPPER_H
#define COST_FUNCTOR_WRAPPER_H

extern "C"
{

void *log_functor_constructor(
	int Mf, int Ma,				// Number of full system and ancilla modes
    const int *full_basis,		// Array of fock states, representing basis of the full system
								// Expected size = fb_size*Mf
	int fb_size,				// Number of fock states in the full basis
	const int *ancilla_basis,	// Array of fock states, representing basis of the ancilla
								// Expected size = ab_size*Ma
	int ab_size,				// Number of fock states in the ancilla basis
	const int *input_state,		// Input fock state. Expected size = Mf
	const int *target_states_focks,	// Two dimensional array of fock states. First index - target state index,
									// second - fock state index in the given target state.
									// Expected size = ts_size1*ts_size2*(Mf-Ma)
	const double *target_states_re, const double *target_states_im,	// Arrays of target states amlitudes (real and imaginary part)
																	// Expected size = ts_size1*ts_size2
	int ts_size1, int ts_size2);	// Number of target states and nonzero amplitudes in target states respectively

void *log_functor_simple_constructor();
void log_functor_destructor(void *functor);
double log_functor_apply(void *functor, const double *x, int x_size);

}

#endif // COST_FUNCTOR_WRAPPER_H

#ifndef COST_FUNCTOR_WRAPPER_H
#define COST_FUNCTOR_WRAPPER_H

extern "C"
{
void *stanisic_functor_constructor(int Mf, int Ma,
                             const int *full_basis,
                             int fb_size,
                             const int *ancilla_basis,
                             int ab_size,
                             const int *input_state,
                             const int *target_states_focks,
                             const double *target_states_re, const double *target_states_im,
                             int ts_size1, int ts_size2);
double stanisic_functor_apply(void *functor, const double *x, int x_size);
}

#endif // COST_FUNCTOR_WRAPPER_H

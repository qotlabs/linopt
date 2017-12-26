#include "cost_functor_wrapper.h"
#include <cost_functor.h>

using namespace linopt;

void *stanisic_functor_constructor(int Mf, int Ma,
                             const int *full_basis,
                             int fb_size,
                             const int *ancilla_basis,
                             int ab_size,
                             const int *input_state,
                             const int *target_states_focks,
                             const double *target_states_re,
                             const double *target_states_im,
                             int ts_size1, int ts_size2)
{
    const int Ms = Mf - Ma;
    basis fb, ab;
    for(auto i = full_basis; i != full_basis + Mf*fb_size; i += Mf)
        fb.insert(fock(i, i + Mf));
    for(auto i = ancilla_basis; i != ancilla_basis + Ma*ab_size; i += Ma)
        ab.insert(fock(i, i + Ma));
    fock is(input_state, input_state + Mf);
    std::vector<state> ts(ts_size1);
    state s;
    fock f(Ms);
    for(int i = 0; i < ts_size1; i++)
    {
        s.clear();
        for(int j = 0; j < ts_size2; j++)
        {
            f = fock(target_states_focks, target_states_focks + Ms);
            target_states_focks += Ms;
            std::complex<double> z(*target_states_re, *target_states_im);
            target_states_re++;
            target_states_im++;
            s.insert(state_element(f, z));
        }
        ts[i] = s;
    }
	stanisic_functor *sf = new stanisic_functor(fb, ab, is, ts);
    return sf;
}

double stanisic_functor_apply(void *functor, const double *x, int x_size)
{
    stanisic_functor *sf = static_cast<stanisic_functor *>(functor);
    Eigen::Map<const point> xmap(x, x_size);
    return sf->operator()(xmap);
}

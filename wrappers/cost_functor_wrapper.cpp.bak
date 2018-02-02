#include "cost_functor_wrapper.h"
#include <vector>
#include <math.h>
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

void *stanisic_functor_simple_constructor()
{
	basis ancilla_basis(2, 4);
	basis full_basis = basis(3, 6)*ancilla_basis;
	fock input_state = {1, 1, 1, 0, 0, 0, 1, 1, 0, 0};
	std::vector<state> B(20);
	B[0][{1, 0, 1, 0, 1, 0}] = M_SQRT1_2;
	B[0][{0, 1, 0, 1, 0, 1}] = M_SQRT1_2;
	B[1][{1, 0, 1, 0, 0, 1}] = M_SQRT1_2;
	B[1][{0, 1, 0, 1, 1, 0}] = M_SQRT1_2;
	B[2][{1, 0, 1, 1, 0, 0}] = M_SQRT1_2;
	B[2][{0, 1, 0, 0, 1, 1}] = M_SQRT1_2;
	B[3][{1, 0, 0, 1, 1, 0}] = M_SQRT1_2;
	B[3][{0, 1, 1, 0, 0, 1}] = M_SQRT1_2;
	B[4][{1, 0, 0, 1, 0, 1}] = M_SQRT1_2;
	B[4][{0, 1, 1, 0, 1, 0}] = M_SQRT1_2;
	B[5][{1, 0, 0, 0, 1, 1}] = M_SQRT1_2;
	B[5][{0, 1, 1, 1, 0, 0}] = M_SQRT1_2;
	B[6][{1, 1, 0, 0, 1, 0}] = M_SQRT1_2;
	B[6][{0, 0, 1, 1, 0, 1}] = M_SQRT1_2;
	B[7][{1, 1, 0, 0, 0, 1}] = M_SQRT1_2;
	B[7][{0, 0, 1, 1, 1, 0}] = M_SQRT1_2;
	B[8][{1, 1, 0, 1, 0, 0}] = M_SQRT1_2;
	B[8][{0, 0, 1, 0, 1, 1}] = M_SQRT1_2;
	B[9][{1, 1, 1, 0, 0, 0}] = M_SQRT1_2;
	B[9][{0, 0, 0, 1, 1, 1}] = M_SQRT1_2;
	B[10][{0, 1, 1, 0, 1, 0}] = M_SQRT1_2;
	B[10][{1, 0, 0, 1, 0, 1}] = M_SQRT1_2;
	B[11][{0, 1, 1, 0, 0, 1}] = M_SQRT1_2;
	B[11][{1, 0, 0, 1, 1, 0}] = M_SQRT1_2;
	B[12][{0, 1, 1, 1, 0, 0}] = M_SQRT1_2;
	B[12][{1, 0, 0, 0, 1, 1}] = M_SQRT1_2;
	B[13][{0, 1, 0, 1, 1, 0}] = M_SQRT1_2;
	B[13][{1, 0, 1, 0, 0, 1}] = M_SQRT1_2;
	B[14][{0, 1, 0, 1, 0, 1}] = M_SQRT1_2;
	B[14][{1, 0, 1, 0, 1, 0}] = M_SQRT1_2;
	B[15][{0, 1, 0, 0, 1, 1}] = M_SQRT1_2;
	B[15][{1, 0, 1, 1, 0, 0}] = M_SQRT1_2;
	B[16][{0, 0, 1, 1, 1, 0}] = M_SQRT1_2;
	B[16][{1, 1, 0, 0, 0, 1}] = M_SQRT1_2;
	B[17][{0, 0, 1, 1, 0, 1}] = M_SQRT1_2;
	B[17][{1, 1, 0, 0, 1, 0}] = M_SQRT1_2;
	B[18][{0, 0, 1, 0, 1, 1}] = M_SQRT1_2;
	B[18][{1, 1, 0, 1, 0, 0}] = M_SQRT1_2;
	B[19][{0, 0, 0, 1, 1, 1}] = M_SQRT1_2;
	B[19][{1, 1, 1, 0, 0, 0}] = M_SQRT1_2;
	stanisic_functor *sf = new stanisic_functor(full_basis, ancilla_basis, input_state, B);
	return sf;
}

void stanisic_functor_destructor(void *functor)
{
    stanisic_functor *sf = static_cast<stanisic_functor *>(functor);
	delete sf;
	return;
}

double stanisic_functor_apply(void *functor, const double *x, int x_size)
{
    stanisic_functor *sf = static_cast<stanisic_functor *>(functor);
    Eigen::Map<const point> xmap(x, x_size);
    return sf->operator()(xmap);
}

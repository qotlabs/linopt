#include "cost_functor.h"

using namespace linopt;

cost_functor::cost_functor(const basis &full_basis,
                           const basis &ancilla_basis,
                           const fock &input_state,
                           const std::vector<state> &target_states):
    target_states(target_states),
    ancilla_basis(ancilla_basis)
{
    C.set_basis(full_basis);
    C.input() = input_state;
}

real_type stanisic_functor::operator()(const unitary_matrix::angles &a)
{
    C.unitary().hurwitz(a);
    state out = C.output_state();
    state postselected;
    real_type p, res = 0.;
    for(auto anc = ancilla_basis.begin(); anc != ancilla_basis.end(); anc++)
    {
        postselected = out.postselect(*anc);
        p = postselected.norm();
        if(p == 0.)
            continue;
        postselected /= p;
        p = p*p;
        for(size_t i = 0; i < target_states.size(); i++)
            res += p * std::pow(std::norm(dot(postselected, target_states[i])), 5);
    }
    return res;
}

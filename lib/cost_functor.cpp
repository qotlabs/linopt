#include "cost_functor.h"

using namespace linopt;

cost_functor::cost_functor(const basis &full_basis,
						   const basis &ancilla_basis,
						   const fock &input_state,
						   const std::vector<state> &target_states):
	target_states(target_states),
	ancilla_basis(ancilla_basis)
{
	C.output_basis() = full_basis;
	C.input_state() = input_state;
}

real_type stanisic_functor::operator()(const point &x)
{
	C.unitary().exp_hermite(x);
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

real_type log_functor::operator()(const point &x)
{
    const real_type epsilon = 1e-2;
    C.unitary().hurwitz(x);
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
		{
            res += epsilon*p/(1. + epsilon - std::norm(dot(postselected, target_states[i])));
		}
	}
	return res;
}

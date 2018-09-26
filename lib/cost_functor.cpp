/* Copyright © 2018, Quantum Optical Technologies Laboratories
 * <https://www.qotlabs.org/en/>
 * Contributed by: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
 *                 Dyakonov Ivan <iv.dyakonov@physics.msu.ru>
 *
 * This file is part of Linopt.
 *
 * Linopt is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Linopt is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Linopt. If not, see <https://www.gnu.org/licenses/>.
 */

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

	C.unitary() = exp_hermite_parametrization(x);
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
	C.unitary() = hurwitz_parametrization(x);
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

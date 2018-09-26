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

#ifndef COST_FUNCTOR_H
#define COST_FUNCTOR_H

#include "circuit.h"
#include <vector>

namespace linopt
{

class cost_functor
{
protected:
	std::vector<state> target_states;
	basis ancilla_basis;
	circuit C;
public:
	cost_functor(const basis &full_basis,
				 const basis &ancilla_basis,
				 const fock &input_state,
				 const std::vector<state> &target_states);
	real_type operator()(const point &x);
};

class stanisic_functor: public cost_functor
{
public:
	using cost_functor::cost_functor;
	real_type operator()(const point &x);
};

class log_functor: public cost_functor
{
public:
	using cost_functor::cost_functor;
	real_type operator()(const point &x);
};

}

#endif // COST_FUNCTOR_H

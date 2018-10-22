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

#include "circuit.h"
#include <functional>

using namespace linopt;

matrix_type &circuit::prepare_uin(const matrix_type &u, const fock &fin)
{
	int tot = fin.total();
	int modes = fin.size();
	int k = 0;
	Uin.resize(u.rows(), tot);
	for(int m = 0; m < modes; m++)
		for(int i = 0; i < fin[m]; i++)
			Uin.col(k++) = u.col(m);
	input_prod_fact = fin.prod_fact();
	return Uin;
}

complex_type circuit::calc_fock_amp(const fock &fout)
{
	int tot = fout.total();
	int modes = fout.size();
	Uinout.resize(tot, Uin.cols());
	int k = 0;
	for(int m = 0; m < modes; m++)
		for(int i = 0; i < fout[m]; i++)
			Uinout.row(k++) = Uin.row(m);
	complex_type perm = permanent(Uinout);
	perm /= std::sqrt(fout.prod_fact() * input_prod_fact);
	return perm;
}

const fock &circuit::get_input_state() const
{
	return _input_state;
}

void circuit::set_input_state(const fock &fin)
{
	output_state_changed = true;
	uin_possibly_changed = true;
	_input_state = fin;
}

const basis circuit::get_output_basis() const
{
	return _output_state.get_basis();
}

void circuit::set_output_basis(const basis &bout)
{
	output_state_changed = true;
	_output_state.set_basis(bout);
}

const matrix_type &circuit::get_unitary() const
{
	return U;
}

void circuit::set_unitary(const matrix_type &u)
{
	output_state_changed = true;
	uin_possibly_changed = true;
	U = u;
}

const state &circuit::output_state()
{
	if(uin_possibly_changed)
	{
		prepare_uin(U, _input_state);
		uin_possibly_changed = false;
	}
	if(output_state_changed)
	{
		_output_state.set_amplitudes(std::bind(&circuit::calc_fock_amp,
									 this, std::placeholders::_1));
		output_state_changed = false;
	}
	return _output_state;
}

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

circuit::circuit():
	uin_possibly_changed(false) {}

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

matrix_type &circuit::unitary()
{
	uin_possibly_changed = true;
	return U;
}

const matrix_type &circuit::unitary() const
{
	return U;
}

fock &circuit::input_state()
{
	uin_possibly_changed = true;
	return _input_state;
}

const fock &circuit::input_state() const
{
	return _input_state;
}

basis &circuit::output_basis()
{
	return _output_basis;
}

const basis &circuit::output_basis() const
{
	return _output_basis;
}

state circuit::output_state()
{
	if(uin_possibly_changed)
	{
		prepare_uin(U, _input_state);
		uin_possibly_changed = false;
	}
	return _output_basis.apply_func(std::bind(&circuit::calc_fock_amp,
											 this, std::placeholders::_1));
}

void circuit::set_input_state(const fock &fin)
{
	this->input_state() = fin;
}

void circuit::set_output_basis(const basis &bout)
{
	this->output_basis() = bout;
}

void circuit::set_unitary(const matrix_type &u)
{
	this->unitary() = u;
}

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

void circuit::copy_columns_on_input(matrix_type &Ucc, const matrix_type &U, const fock &fin)
{
	const int tot = fin.total();
	const int modes = fin.size();
	int k = 0;
	Ucc.resize(U.rows(), tot);
	for(int m = 0; m < modes; m++)
		for(int i = 0; i < fin[m]; i++)
			Ucc.col(k++) = U.col(m);
}

void circuit::copy_rows_on_output(matrix_type &Ucr, const matrix_type &U, const fock &fout)
{
	const int tot = fout.total();
	const int modes = fout.size();
	Ucr.resize(tot, U.cols());
	int k = 0;
	for(int m = 0; m < modes; m++)
		for(int i = 0; i < fout[m]; i++)
			Ucr.row(k++) = U.row(m);
}

complex_type circuit::calc_fock_amp(const fock &fout) const
{
	matrix_type Uout;
	copy_rows_on_output(Uout, unitary, fout);
	const auto tot = fout.total();
	const auto out_prod_fact = fout.prod_fact();
	matrix_type Uoutin(tot, tot);
	complex_type amp = 0.;
	for(const auto &in_elem: input_state)
	{
		const fock &fin = in_elem.first;
		if(fin.total() != tot)
			continue;
		copy_columns_on_input(Uoutin, Uout, fin);
		amp += permanent(Uoutin) * in_elem.second /
								std::sqrt(out_prod_fact * fin.prod_fact());
	}
	return amp;
}

complex_type circuit::calc_fock_amp_1(const uin_fin &precomputed, const fock &fout) const
{
	if(fout.total()	!= precomputed.tot)
		return 0.;
	matrix_type Uinout;
	copy_rows_on_output(Uinout, precomputed.Uin, fout);
	return permanent(Uinout) * precomputed.mult / std::sqrt(fout.prod_fact());
}

const state &circuit::get_input_state() const
{
	return input_state;
}

void circuit::set_input_state(const state &s)
{
	output_state_valid = false;
	input_state = s;
}

const basis circuit::get_output_basis() const
{
	return _output_state.get_basis();
}

void circuit::set_output_basis(const basis &bout)
{
	output_state_valid = false;
	_output_state.set_basis(bout);
}

const matrix_type &circuit::get_unitary() const
{
	return unitary;
}

void circuit::set_unitary(const matrix_type &U)
{
	output_state_valid = false;
	unitary = U;
}

template<class exec_policy>
const state &circuit::output_state()
{
	using namespace std;
	using namespace std::placeholders;
	if(!output_state_valid)
	{
		if(input_state.size() > 1)
		{
			_output_state.set_amplitudes<exec_policy>(
				bind(&circuit::calc_fock_amp, this, _1));
		}
		else
		{
			uin_fin precomputed;
			const fock &fin = input_state.begin()->first;
			const auto &amp = input_state.begin()->second;
			copy_columns_on_input(precomputed.Uin, unitary, fin);
			precomputed.tot = fin.total();
			precomputed.mult = amp / std::sqrt(fin.prod_fact());
			_output_state.set_amplitudes<exec_policy>(
				bind(&circuit::calc_fock_amp_1, this, ref(precomputed), _1));
		}
		output_state_valid = true;
	}
	return _output_state;
}

template const state &circuit::output_state<execution::seq>();
template const state &circuit::output_state<execution::par>();

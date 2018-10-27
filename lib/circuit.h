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

#ifndef CIRCUIT_H
#define CIRCUIT_H

#include "states.h"
#include "matrix.h"

namespace linopt
{

class circuit
{
private:
	state input_state;
	matrix_type unitary;
	state _output_state;

	// Cache

	// Indicates whether variable `_output_state` is valid
	bool output_state_valid = true;

	// Cache struct for `calc_fock_amp_1()`
	struct uin_fin
	{
		matrix_type Uin;
		int tot;
		complex_type mult;
	};

	// Members

	static void copy_columns_on_input(matrix_type &Ucc, const matrix_type &U, const fock &fin);
	static void copy_rows_on_output(matrix_type &Ucr, const matrix_type &Uin, const fock &fout);
	// Calculates amplitude corresponding to the output Fock state `fout`
	complex_type calc_fock_amp(const fock &fout) const;
	// Optimized version of `calc_fock_amp()` when `input_state` has size = 1
	complex_type calc_fock_amp_1(const uin_fin &precomputed, const fock &fout) const;

public:
	const state &get_input_state() const;
	void set_input_state(const state &s);
	const basis get_output_basis() const;
	void set_output_basis(const basis &bout);
	const matrix_type &get_unitary() const;
	void set_unitary(const matrix_type &U);
	template<class exec_policity = execution::seq>
	const state &output_state();
};

}

#endif // CIRCUIT_H

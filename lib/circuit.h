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
	matrix_type U;
	matrix_type Uin;
	bool uin_possibly_changed;
	matrix_type Uinout;
	fock _input_state;
	real_type input_prod_fact;
	basis _output_basis;
	matrix_type &prepare_uin(const matrix_type &u, const fock &fin);
	complex_type calc_fock_amp(const fock &fout);

public:
	circuit();
	matrix_type &unitary();
	const matrix_type &unitary() const;
	fock &input_state();
	const fock &input_state() const;
	basis &output_basis();
	const basis &output_basis() const;
	state output_state();
	void set_input_state(const fock &fin);
	void set_output_basis(const basis &bout);
	void set_unitary(const matrix_type &u);
};

}

#endif // CIRCUIT_H

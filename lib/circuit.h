/* Copyright © 2018-2020, Quantum Optical Technologies Laboratories
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

#ifndef _LINOPT_CIRCUIT_H
#define _LINOPT_CIRCUIT_H

#include "states.h"
#include "matrix.h"

namespace linopt
{

class Circuit
{
public:
	const State &getInputState() const;
	void setInputState(const State &s);
	const Basis getOutputBasis() const;
	void setOutputBasis(const Basis &bout);
	const Matrix &getUnitary() const;
	void setUnitary(const Matrix &U);
	template<typename ExecPolicy = execution::Seq>
	const State &outputState();

private:
	State inputState;
	Matrix unitary;
	State outputState_;

	// Indicates whether variable `outputState_` is valid
	bool outputStateValid = true;

	// Cache struct for `calcFockAmp1()`
	struct UinFin
	{
		Matrix Uin;
		int tot;
		Complex mult;
	};

	static void copyColumnsOnInput(Matrix &Ucc, const Matrix &U, const Fock &fin);
	static void copyRowsOnOutput(Matrix &Ucr, const Matrix &Uin, const Fock &fout);
	// Calculates amplitude corresponding to the output Fock state `fout`
	Complex calcFockAmp(const Fock &fout) const;
	// Optimized version of `calcFockAmp()` when `inputState` has size = 1
	Complex calcFockAmp1(const UinFin &precomputed, const Fock &fout) const;
};

} // Namespace linopt

#endif // _LINOPT_CIRCUIT_H

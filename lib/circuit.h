// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

#pragma once

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
	template<typename ExecPolicy = std::execution::sequenced_policy>
	const State &outputState();

private:
	State inputState;
	Matrix unitary;
	State outputState_;

	/// Whether variable `outputState_` is valid
	bool outputStateValid = true;

	/// Cache struct for `calcFockAmp1()`
	struct UinFin
	{
		Matrix Uin;
		int tot;
		Complex mult;
	};

	static void copyColumnsOnInput(Matrix &Ucc, const Matrix &U, const Fock &fin);
	static void copyRowsOnOutput(Matrix &Ucr, const Matrix &Uin, const Fock &fout);

	/// Calculates amplitude corresponding to the output Fock state `fout`
	Complex calcFockAmp(const Fock &fout) const;

	/// Optimized version of `calcFockAmp()` when `inputState` has size = 1
	Complex calcFockAmp1(const UinFin &precomputed, const Fock &fout) const;
};

} // namespace linopt

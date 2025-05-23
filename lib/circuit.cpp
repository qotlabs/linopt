// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

#include "circuit.h"
#include <functional>

using namespace linopt;

void Circuit::copyColumnsOnInput(Matrix &Ucc, const Matrix &U, const Fock &fin)
{
	const int tot = fin.total();
	const int modes = fin.size();
	int k = 0;
	Ucc.resize(U.rows(), tot);
	for(int m = 0; m < modes; m++)
		for(int i = 0; i < fin[m]; i++)
			Ucc.col(k++) = U.col(m);
}

void Circuit::copyRowsOnOutput(Matrix &Ucr, const Matrix &U, const Fock &fout)
{
	const int tot = fout.total();
	const int modes = fout.size();
	Ucr.resize(tot, U.cols());
	int k = 0;
	for(int m = 0; m < modes; m++)
		for(int i = 0; i < fout[m]; i++)
			Ucr.row(k++) = U.row(m);
}

Complex Circuit::calcFockAmp(const Fock &fout) const
{
	Matrix Uout;
	copyRowsOnOutput(Uout, unitary, fout);
	const auto tot = fout.total();
	const auto outProdFact = fout.prodFact();
	Matrix Uoutin(tot, tot);
	Complex amp = 0.;
	for(const auto &inElem: inputState)
	{
		const Fock &fin = inElem.first;
		if(fin.total() != tot)
			continue;
		copyColumnsOnInput(Uoutin, Uout, fin);
		amp += permanent(Uoutin) * inElem.second /
								std::sqrt(outProdFact * fin.prodFact());
	}
	return amp;
}

Complex Circuit::calcFockAmp1(const UinFin &precomputed, const Fock &fout) const
{
	if(fout.total()	!= precomputed.tot)
		return 0.;
	Matrix Uinout;
	copyRowsOnOutput(Uinout, precomputed.Uin, fout);
	return permanent(Uinout) * precomputed.mult / std::sqrt(fout.prodFact());
}

const State &Circuit::getInputState() const
{
	return inputState;
}

void Circuit::setInputState(const State &s)
{
	outputStateValid = false;
	inputState = s;
}

const Basis Circuit::getOutputBasis() const
{
	return outputState_.getBasis();
}

void Circuit::setOutputBasis(const Basis &bout)
{
	outputStateValid = false;
	outputState_.setBasis(bout);
}

const Matrix &Circuit::getUnitary() const
{
	return unitary;
}

void Circuit::setUnitary(const Matrix &U)
{
	outputStateValid = false;
	unitary = U;
}

template<typename ExecPolicy>
const State &Circuit::outputState()
{
	using namespace std;
	using namespace std::placeholders;
	if(!outputStateValid)
	{
		if(inputState.size() > 1)
		{
			outputState_.setAmplitudes<ExecPolicy>(
				bind(&Circuit::calcFockAmp, this, _1));
		}
		else
		{
			UinFin precomputed;
			const Fock &fin = inputState.begin()->first;
			const auto &amp = inputState.begin()->second;
			copyColumnsOnInput(precomputed.Uin, unitary, fin);
			precomputed.tot = fin.total();
			precomputed.mult = amp / std::sqrt(fin.prodFact());
			outputState_.setAmplitudes<ExecPolicy>(
				bind(&Circuit::calcFockAmp1, this, ref(precomputed), _1));
		}
		outputStateValid = true;
	}
	return outputState_;
}

template const State &Circuit::outputState<std::execution::sequenced_policy>();
template const State &Circuit::outputState<std::execution::parallel_policy>();

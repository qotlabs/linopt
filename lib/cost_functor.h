// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

#pragma once

#include "circuit.h"
#include <vector>

namespace linopt
{

class CostFunctor
{
public:
	CostFunctor(const Basis &fullBasis,
				const Basis &ancillaBasis,
				const Fock &inputState,
				const std::vector<State> &targetStates);
	Real operator()(const Point &x);

protected:
	std::vector<State> targetStates;
	Basis ancillaBasis;
	Circuit C;
};

class StanisicFunctor: public CostFunctor
{
public:
	using CostFunctor::CostFunctor;
	Real operator()(const Point &x);
};

class LogFunctor: public CostFunctor
{
public:
	using CostFunctor::CostFunctor;
	Real operator()(const Point &x);
};

} // namespace linopt

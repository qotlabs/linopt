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

#ifndef _LINOPT_COST_FUNCTOR_H
#define _LINOPT_COST_FUNCTOR_H

#include <vector>
#include "circuit.h"

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

} // Namespace linopt

#endif // _LINOPT_COST_FUNCTOR_H

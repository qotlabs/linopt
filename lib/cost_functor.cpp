/* Copyright © 2018, 2019, Quantum Optical Technologies Laboratories
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

#include "cost_functor.h"

using namespace linopt;

CostFunctor::CostFunctor(const Basis &fullBasis,
						 const Basis &ancillaBasis,
						 const Fock &inputState,
						 const std::vector<State> &targetStates):
	targetStates(targetStates),
	ancillaBasis(ancillaBasis)
{
	C.setOutputBasis(fullBasis);
	C.setInputState(inputState);
}

Real StanisicFunctor::operator()(const Point &x)
{
	C.setUnitary(expHermiteParametrization(x));
	State out = C.outputState();
	State postselected;
	Real p, res = 0.;
	for(auto anc = ancillaBasis.begin(); anc != ancillaBasis.end(); anc++)
	{
		postselected = out.postselect(*anc);
		p = postselected.norm();
		if(p == 0.)
			continue;
		postselected /= p;
		p = p*p;
		for(size_t i = 0; i < targetStates.size(); i++)
			res += p * std::pow(std::norm(dot(postselected, targetStates[i])), 5);
	}
	return res;
}

Real LogFunctor::operator()(const Point &x)
{
	const Real epsilon = 1e-2;
	C.setUnitary(hurwitzParametrization(x));
	State out = C.outputState();
	State postselected;
	Real p, res = 0.;
	for(auto anc = ancillaBasis.begin(); anc != ancillaBasis.end(); anc++)
	{
		postselected = out.postselect(*anc);
		p = postselected.norm();
		if(p == 0.)
			continue;
		postselected /= p;
		p = p*p;
		for(size_t i = 0; i < targetStates.size(); i++)
		{
			res += epsilon*p/(1. + epsilon - std::norm(dot(postselected, targetStates[i])));
		}
	}
	return res;
}

// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

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

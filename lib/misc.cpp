// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

#include "misc.h"

using namespace linopt;

std::ostream& linopt::printComplex(std::ostream &stream, const Complex &x)
{
	const Real re = x.real();
	const Real im = x.imag();
	stream << re;
	if(im >= 0)
		stream << "+";
	stream << im << "j";
	return stream;
}

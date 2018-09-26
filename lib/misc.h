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

#ifndef MISC_H
#define MISC_H

#include <cmath>
#include <ostream>
#include "types.h"

namespace linopt
{

template<typename T>
inline T mod(T x, T a)
{
	x = std::fmod(x, a);
	return (x < 0) ? x + a : x;
}

template<typename T>
std::ostream& print_array(std::ostream &stream, const T &a,
						  const char *b1 = "{",
						  const char *delim = ", ",
						  const char *b2 = "}")
{
	if(a.empty())
	{
		stream << b1 << b2;
	}
	else
	{
		stream << b1;
		for(auto iter = a.begin(); iter != --a.end(); iter++)
			stream << *iter << delim;
		stream << *(--a.end()) << b2;
	}
	return stream;
}

std::ostream& print_complex(std::ostream &stream, const complex_type &x);

}

#endif // MISC_H

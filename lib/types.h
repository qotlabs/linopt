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

#ifndef _LINOPT_TYPES_H
#define _LINOPT_TYPES_H

#include <complex>

namespace linopt
{

using Real = double;
using Complex = std::complex<Real>;

/** @ingroup matrix
 * @brief Default precision for numeric comparison operations.
 */
constexpr Real defaultEpsilon = 1e-15;

/**
 * @brief The namespace containing execution policies used by several library
 * functions.
 */
namespace execution
{
	class Seq {};	///< Sequential execution policy.
	class Par {};	///< Parallel execution policy.
} // Namespace execution

} // Namespace linopt

#endif // _LINOPT_TYPES_H

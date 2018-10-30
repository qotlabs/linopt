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

/** @mainpage
 * Linopt is a library for linear optics calculations. It consists of C++11 core
 * library and Python bindings to it.
 */

#ifndef LINOPT_H
#define LINOPT_H

#include "types.h"
#include "matrix.h"
#include "states.h"
#include "circuit.h"
#include "circuit_design.h"

/**
 * @brief The main namespace containing all library classes, functions, etc.
 */
namespace linopt
{

/**
  * @brief Library version specification.
  */
constexpr struct version
{
	static constexpr auto string = "0.2.0";
	static constexpr int major = 0;
	static constexpr int minor = 2;
	static constexpr int patch = 0;
} version;

}

#endif // LINOPT_H

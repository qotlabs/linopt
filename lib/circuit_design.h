/* Copyright © 2018, Quantum Optical Technologies Laboratories
 * <https://www.qotlabs.org/en/>
 * Contributed by: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
 *                 Dyakonov Ivan <iv.dyakonov@physics.msu.ru>
 *				   Fldjyan Suren <fldzhian.sa17@physics.msu.ru>
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

#ifndef CIRCUIT_DESIGN_H
#define CIRCUIT_DESIGN_H

#include "linopt.h"
#include "matrix.h"

namespace linopt
{

void clements_design(matrix_type &M, const point &x, const point &y);
void clements_design(matrix_type &M, const point &x);
matrix_type clements_design(const point &x, const point &y);
matrix_type clements_design(const point &x);
void get_clements_design(matrix_type &M, point &x, real_type eps  = default_epsilon);
point get_clements_design(matrix_type &M, real_type eps = default_epsilon);

}

#endif // CIRCUIT_DESIGN_H

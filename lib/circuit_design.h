/* Copyright © 2018, 2019, Quantum Optical Technologies Laboratories
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

#include "matrix.h"

namespace linopt
{

void clementsDesign(Matrix &M, const Point &x, const Point &y);
void clementsDesign(Matrix &M, const Point &x);
Matrix clementsDesign(const Point &x, const Point &y);
Matrix clementsDesign(const Point &x);
void getClementsDesign(Matrix &M, Point &x, Real eps  = defaultEpsilon);
Point getClementsDesign(Matrix &M, Real eps = defaultEpsilon);

} // Namespace linopt

#endif // CIRCUIT_DESIGN_H

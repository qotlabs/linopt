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

#ifndef MATRIX_H
#define MATRIX_H

#include <Eigen/Dense>
#include <vector>
#include "types.h"

namespace linopt
{

using Matrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using Column = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;
using Row = Eigen::Matrix<Complex, 1, Eigen::Dynamic>;
using Point = std::vector<Real>;

void hurwitzParametrization(Matrix &M, const Point &x);
Matrix hurwitzParametrization(const Point &x);
void expHermiteParametrization(Matrix &M, const Point &x);
Matrix expHermiteParametrization(const Point &x);

Real matrixFidelity(const Matrix &A, const Matrix &B);

bool isColumnUnitary(const Matrix &M, Real eps = defaultEpsilon);
bool isRowUnitary(const Matrix &M, Real eps = defaultEpsilon);
bool isUnitary(const Matrix &M, Real eps = defaultEpsilon);

Complex permanent(const Matrix &M);

} // Namespace linopt

#endif // MATRIX_H

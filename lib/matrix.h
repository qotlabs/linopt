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

#ifndef MATRIX_H
#define MATRIX_H

#include <Eigen/Dense>
#include <vector>
#include "linopt.h"
#include "types.h"

namespace linopt
{

typedef Eigen::Matrix<complex_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> matrix_type;
typedef Eigen::Matrix<complex_type, Eigen::Dynamic, 1> vector_type;
typedef std::vector<real_type> point;

/** @ingroup matrix
 * @brief Default precision for numeric comparison operations
 */
const real_type default_epsilon = 1e-15;

void hurwitz_parametrization(matrix_type &M, const point &x);
matrix_type hurwitz_parametrization(const point &x);
void exp_hermite_parametrization(matrix_type &M, const point &x);
matrix_type exp_hermite_parametrization(const point &x);

real_type matrix_fidelity(const matrix_type &A, const matrix_type &B);

bool is_column_unitary(const matrix_type &M, real_type eps = default_epsilon);
bool is_row_unitary(const matrix_type &M, real_type eps = default_epsilon);
bool is_unitary(const matrix_type &M, real_type eps = default_epsilon);

complex_type permanent(const matrix_type &M);
complex_type ppermanent(const matrix_type &M);	// TODO: Parallel version

}

#endif // MATRIX_H

#ifndef MATRIX_H
#define MATRIX_H

#include <Eigen/Dense>
#include "types.h"
#include "states.h"

namespace linopt
{

const real_type default_epsilon = 1e-15;

void hurwitz_parametrization(matrix_type &M, const point &x);
matrix_type hurwitz_parametrization(const point &x);
void exp_hermite_parametrization(matrix_type &M, const point &x);
matrix_type exp_hermite_parametrization(const point &x);

bool is_column_unitary(const matrix_type &M, real_type eps = default_epsilon);
bool is_row_unitary(const matrix_type &M, real_type eps = default_epsilon);
bool is_unitary(const matrix_type &M, real_type eps = default_epsilon);

complex_type permanent(const matrix_type &M);
complex_type ppermanent(const matrix_type &M);	// TODO: Parallel version

}

#endif // MATRIX_H

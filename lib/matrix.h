#ifndef MATRIX_H
#define MATRIX_H

#include <Eigen/Dense>
#include "states.h"

namespace linopt
{

typedef Eigen::Matrix<complex_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> matrix_type;
typedef Eigen::Matrix<complex_type, Eigen::Dynamic, 1> vector_type;
typedef std::vector<real_type> point;

const real_type default_epsilon = 1e-15;

void hurwitz_parametrization(matrix_type &M, const point &x);
matrix_type hurwitz_parametrization(const point &x);
void exp_hermite_parametrization(matrix_type &M, const point &x);
matrix_type exp_hermite_parametrization(const point &x);

bool is_column_unitary(const matrix_type &M, real_type eps = default_epsilon);
bool is_row_unitary(const matrix_type &M, real_type eps = default_epsilon);
bool is_unitary(const matrix_type &M, real_type eps = default_epsilon);

complex_type permanent(const matrix_type &M);

template<typename T>
static inline T mod(T x, T a)
{
	x = std::fmod(x, a);
	return (x < 0) ? x + a : x;
}

static inline real_type pyramid(real_type x, real_type a)
{
	x = std::fabs(x);
	x = ((unsigned)(x/a) & 1) ? -x : x;
	return mod(x, a);
}

}

#endif // MATRIX_H

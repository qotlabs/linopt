#ifndef MATRIX_H
#define MATRIX_H

#include <eigen3/Eigen/Dense>
#include "states.h"

namespace linopt
{

typedef Eigen::Matrix<complex_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrix;

complex_type permanent(const matrix &M);

}

#endif // MATRIX_H

#ifndef MATRIX_H
#define MATRIX_H

#include <eigen3/Eigen/Dense>
#include "states.h"

namespace linopt
{

typedef Eigen::Matrix<complex_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> matrix_type;
typedef Eigen::Matrix<complex_type, Eigen::Dynamic, 1> vector_type;

class unitary_matrix: public matrix_type
{
public:
    typedef std::vector<real_type> angles;
    using matrix_type::matrix_type;
    unitary_matrix &hurwitz(const angles &a);
};

complex_type permanent(const matrix_type &M);

}

#endif // MATRIX_H

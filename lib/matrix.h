#ifndef MATRIX_H
#define MATRIX_H

#include <Eigen/Dense>
#include "states.h"

namespace linopt
{

typedef Eigen::Matrix<complex_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> matrix_type;
typedef Eigen::Matrix<complex_type, Eigen::Dynamic, 1> vector_type;
typedef Eigen::Matrix<real_type, Eigen::Dynamic, 1> point;

class unitary_matrix: public matrix_type
{
public:
    static constexpr real_type default_epsilon = 1e-15;
    using matrix_type::matrix_type;
    unitary_matrix &hurwitz(const point &x);
    unitary_matrix &exp_hermite(const point &x);
    bool is_column_unitary(real_type eps = default_epsilon) const;
    bool is_row_unitary(real_type eps = default_epsilon) const;
    bool is_unitary(real_type eps = default_epsilon) const;
};

complex_type permanent(const matrix_type &M);

static inline real_type mod(real_type x, real_type a)
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

#ifndef TYPES_H
#define TYPES_H

#include <functional>
#include <complex>
#include <vector>
#include <Eigen/Dense>

namespace linopt
{

typedef double real_type;
typedef std::complex<real_type> complex_type;

class fock;
class basis;
class state;

typedef std::pair<fock, complex_type> state_element;
typedef std::function<complex_type(const fock&)> basis_func;

typedef Eigen::Matrix<complex_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> matrix_type;
typedef Eigen::Matrix<complex_type, Eigen::Dynamic, 1> vector_type;
typedef std::vector<real_type> point;

}

#endif // TYPES_H

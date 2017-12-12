#include "matrix.h"

using namespace linopt;

complex_type linopt::permanent(const matrix &M)
{
    int nmax = 1 << (M.rows() - 1);
    real_type mult, sign = 1;
    complex_type perm = 0.;
    Eigen::VectorXcd sum = M.colwise().sum();
    for(int n = 0; n < nmax; n++)
    {
        perm += sign * sum.prod();
        sign *= -1;
        int i = n ^ (n + 1);
        i++;
        mult = ((n + 1) & i) ? 2 : -2;
        i = __builtin_ctz(i) - 1;
        sum += mult * M.row(i);
    }
    perm /= nmax;
    return perm;
}

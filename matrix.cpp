#include <cassert>
#include <limits.h>
#include "matrix.h"

using namespace linopt;

complex_type linopt::permanent(const matrix &M)
{
    assert(M.rows() == M.cols());
    if(M.rows() == 2)
        return M(0, 0)*M(1, 1) + M(0, 1)*M(1, 0);
    if(M.rows() == 3)
        return M(0, 0)*(M(1, 1)*M(2, 2) + M(1, 2)*M(2, 1)) +
               M(0, 1)*(M(1, 0)*M(2, 2) + M(1, 2)*M(2, 0)) +
               M(0, 2)*(M(1, 0)*M(2, 1) + M(1, 1)*M(2, 0));
    if(M.rows() == 4)
        return M(0,1)*M(1,3)*M(2,2)*M(3,0) + M(0,1)*M(1,2)*M(2,3)*M(3,0) + M(0,0)*M(1,3)*M(2,2)*M(3,1) +
               M(0,0)*M(1,2)*M(2,3)*M(3,1) + M(0,1)*M(1,3)*M(2,0)*M(3,2) + M(0,0)*M(1,3)*M(2,1)*M(3,2) +
               M(0,1)*M(1,0)*M(2,3)*M(3,2) + M(0,0)*M(1,1)*M(2,3)*M(3,2) + M(0,3)*(M(1,2)*(M(2,1)*M(3,0) +
               M(2,0)*M(3,1)) + M(1,1)*(M(2,2)*M(3,0) + M(2,0)*M(3,2)) + M(1,0)*(M(2,2)*M(3,1) +
               M(2,1)*M(3,2))) + M(0,1)*M(1,2)*M(2,0)*M(3,3) + M(0,0)*M(1,2)*M(2,1)*M(3,3) +
               M(0,1)*M(1,0)*M(2,2)*M(3,3) + M(0,0)*M(1,1)*M(2,2)*M(3,3) + M(0,2)*(M(1,3)*(M(2,1)*M(3,0) +
               M(2,0)*M(3,1)) + M(1,1)*(M(2,3)*M(3,0) + M(2,0)*M(3,3)) + M(1,0)*(M(2,3)*M(3,1) +
               M(2,1)*M(3,3)));
    assert(static_cast<unsigned long>(M.rows()) <= CHAR_BIT*sizeof(unsigned long));
    complex_type perm = 0.;
    Eigen::VectorXcd sum = M.colwise().sum();
    real_type mult;
    unsigned long n = 0ul;
    unsigned long nmax = 1ul << (M.rows() - 1ul);
    while(n < nmax)
    {
        complex_type sp = sum.prod();
        perm += (n & 1ul) ? -sp : sp;
        unsigned long i = n ^ (n + 1ul);
        i++;
        n++;
        mult = (n & i) ? 2 : -2;
        i = __builtin_ctz(i) - 1ul;
        sum += mult * M.row(i);
    }
    perm /= nmax;
    return perm;
}

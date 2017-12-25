#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <cstdlib>
#include <math.h>
#include <ctime>

#include "linopt.h"
#include "bfgs.h"
#include "hco.h"

using namespace std;
using namespace linopt;

int srandomdev(void)
{
    int fd, seed;
    fd = open("/dev/urandom", O_RDONLY);
    if(read(fd, &seed, sizeof(int)) != -1)
        srand(seed);
    else
        srand((unsigned)time(0));
    close(fd);
    return seed;
}

int main()
{
    srandomdev();
    basis full_basis = basis(2, 4)*basis(2, 4);
    basis ancilla_basis = basis(2, 4);
    fock input_state = {1, 1, 1, 1, 0, 0, 0, 0};
    vector<state> B(6);
    B[0][{1, 0, 1, 0}] =  M_SQRT1_2;
    B[0][{0, 1, 0, 1}] =  M_SQRT1_2;
    B[1][{1, 0, 1, 0}] =  M_SQRT1_2;
    B[1][{0, 1, 0, 1}] = -M_SQRT1_2;
    B[2][{1, 0, 0, 1}] =  M_SQRT1_2;
    B[2][{0, 1, 1, 0}] =  M_SQRT1_2;
    B[3][{1, 0, 0, 1}] =  M_SQRT1_2;
    B[3][{0, 1, 1, 0}] = -M_SQRT1_2;
    B[4][{1, 1, 0, 0}] =  M_SQRT1_2;
    B[4][{0, 0, 1, 1}] =  M_SQRT1_2;
    B[5][{1, 1, 0, 0}] =  M_SQRT1_2;
    B[5][{0, 0, 1, 1}] = -M_SQRT1_2;
    point a(64);
    a = 0.5 * a.setRandom().array() + 0.5;
    stanisic_functor cf(full_basis, ancilla_basis, input_state, B);
    stop_criterion crit(2000, 1e-4, 0, 1e-8, 0);
    real_type val = bfgs(cf, a, crit);
    cout << val << endl;
    return 0;
}

#include <iostream>
#include <unistd.h>
#include <fcntl.h>

#include "linopt.h"
#include "optimization.h"
#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"
#include <cppoptlib/solver/bfgssolver.h>
#include <cppoptlib/solver/neldermeadsolver.h>
#include <cppoptlib/solver/conjugatedgradientdescentsolver.h>

using namespace std;
using namespace linopt;
using namespace cppoptlib;

int srandomdev(void)
{
    int fd, seed;
    fd = open("/dev/urandom", O_RDONLY);
    if(read(fd, &seed, sizeof(int)) != -1)
        srandom(seed);
    else
        srandom((unsigned)time(0));
    close(fd);
    return seed;
}

class problem : public Problem<real_type>
{
private:
    log_functor f;
public:
    problem(const basis &full_basis, const basis &ancilla_basis,
                     const fock &input_state, const vector<state> &target_states):
        f(full_basis, ancilla_basis, input_state, target_states) {};
    real_type value(const TVector &x)
    {
        real_type fval = f(x);
        cout << fval << endl;
        return -fval;
    }
};

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
    stop_criterion crit(1000, 1e-2, 0, 1e-8, 0);
    real_type val = bfgs(cf, a, crit);
    cout << val << endl;
    return 0;
}

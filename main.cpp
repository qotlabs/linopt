#include <iostream>
#include <unistd.h>
#include <fcntl.h>

#include "linopt.h"
#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"
#include <cppoptlib/solver/bfgssolver.h>
#include <cppoptlib/solver/neldermeadsolver.h>

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

class stanisic_problem : public Problem<real_type>
{
private:
    stanisic_functor sf;
public:
    stanisic_problem(const basis &full_basis, const basis &ancilla_basis,
                     const fock &input_state, const vector<state> &target_states):
        sf(full_basis, ancilla_basis, input_state, target_states) {};
    real_type value(const TVector &x)
    {
        real_type f = sf(x);
        return -f;
    }
};

int main()
{
    srandomdev();
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
    stanisic_problem cost_function(basis(2, 4)*basis(2, 4), basis(2, 4), input_state, B);
    BfgsSolver<stanisic_problem> solver;
    point a(64);
    a = 0.5 * a.setRandom().array() + 0.5;
    solver.minimize(cost_function, a);
    cout << -cost_function(a) << endl;
    return 0;
}

#include <iostream>

#include "linopt.h"
#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"
#include <cppoptlib/solver/bfgssolver.h>
#include <cppoptlib/solver/neldermeadsolver.h>

using namespace std;
using namespace linopt;
using namespace cppoptlib;

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
        unitary_matrix::angles a(x.data(), x.data() + x.size());
        real_type f = sf(a);
        //cout << f << endl;
        return -f;
    }
};

int main()
{
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
    srand((unsigned int) time(0));
    stanisic_problem cost_function(basis(2, 4)*basis(2, 4), basis(2, 4), input_state, B);
    BfgsSolver<stanisic_problem> solver;
    Eigen::VectorXd a(64);
    a.setRandom();
    solver.minimize(cost_function, a);
    cout << -cost_function(a) << endl;
    return 0;
}

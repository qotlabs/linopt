#include <iostream>

#include "linopt.h"

using namespace std;
using namespace linopt;

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
    unitary_matrix::angles a(64, 0.4);
    basis b(4, 8);
    chip C;
    C.unitary().hurwitz(a);
    cout << C.unitary() << endl;
    C.input() = input_state;
    C.set_basis(b);
    cout << b.size() << endl;
    matrix_type m(4, 4);
    m.setRandom();
    complex_type z = 0;
    for(int i = 0; i < 3300000; i++)
    {
        z += permanent(m);
        //C.output_state();
    }
    cout << C.output_state() << endl;
    return 0;
}

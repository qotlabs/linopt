#include <iostream>

#include "linopt.h"

using namespace std;
using namespace linopt;

int main()
{
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
    matrix m = matrix::Random(4, 4);
    for(int i = 0; i < 10000000; i++)
    {
        permanent(m);
    }
    cout << permanent(m) << endl;
    return 0;
}

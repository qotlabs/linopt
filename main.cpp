#include <iostream>

#include "linopt.h"

using namespace std;
using namespace linopt;

complex_type f(fock x)
{
    return complex_type(x[0]*x[0], x[1]*x[1]);
}

int main()
{
    basis b1 = {{1, 2, 1, 5}, {2, 3, 1, 5}};
    basis b2 = {{3, 4, 8, 5}, {4, 3, 1, 5}};
    basis b3 = b1 + b2;
    cout << b1 << endl;
    cout << b2 << endl;
    cout << b3 << endl;
    state s = b3.apply_func(f);
    cout << s << endl;
    fock f1 = {1, 2};
    fock f2 = {1 ,3};
    cout << (f1 < f2) << endl;
    cout << b3.postselect({3, 1, 5}) << endl;
    cout << s.postselect({3, 1, 5}) << endl;
    cout << s.get_basis() << endl;
    cout << s.norm() << endl;
    cout << s.normalize() << endl;
    cout << s << endl;
    return 0;
}

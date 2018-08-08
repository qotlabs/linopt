#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <cstdlib>
#include <math.h>
#include <ctime>

#include <linopt.h>

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

real_type ackley(const point &x)
{
    real_type d = 0.;
    real_type s = 0.;
    for(int i = 0; i < x.size(); i++)
    {
        d += x[i]*x[i];
        s += std::cos(2.*M_PI*x[i]);
    }
    s /= x.size();
    d /= x.size();
    return 20.*std::exp(-0.2*std::sqrt(d)) + std::exp(s) - 20 - M_E;
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
    //a = 0.5 * a.setRandom().array() + 0.5;
	a = 10. * a.setRandom();
	chip C;
    C.unitary().hurwitz(a);
    C.input_state() = input_state;
	C.output_basis() = full_basis;
	state out = C.output_state();
    for(auto anc = ancilla_basis.begin(); anc != ancilla_basis.end(); anc++)
		cout << out.postselect(*anc).normalize() << endl;
    return 0;
}

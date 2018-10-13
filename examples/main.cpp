#include <iostream>
#include <random>

#include <linopt.h>

using namespace std;
using namespace linopt;

int main()
{
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
	random_device rand_dev;
	mt19937 rand_gen(rand_dev());
	uniform_real_distribution<real_type> distribution(0., 1.);
	for(unsigned i = 0; i < a.size(); i++)
		a[i] = distribution(rand_gen);
	circuit C;
	C.set_unitary(hurwitz_parametrization(a));
	C.set_input_state(input_state);
	C.set_output_basis(full_basis);
	state out = C.output_state();
	for(const auto &anc : ancilla_basis)
		cout << out.postselect(anc).normalize() << endl;
	return EXIT_SUCCESS;
}

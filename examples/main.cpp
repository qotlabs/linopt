#include <iostream>
#include <random>

#include <states.h>
#include <circuit.h>

using namespace std;
using namespace linopt;

int main()
{
	const int nphot = 4;
	const int modes = 8;
	basis full_basis = basis(nphot/2, modes/2)*basis(nphot/2, modes/2);
	fock input_state(modes);
	for(int i = 0; i < modes; i++)
		input_state[i] = 1;
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
	point a(modes*modes);
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

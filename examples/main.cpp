#include <iostream>
#include <random>

#include <linopt.h>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace linopt;

int main()
{
	const int nphot = 10;
	const int modes = 20;
	basis full_basis = basis(nphot/2, modes/2)*basis(nphot/2, modes/2);
	basis ancilla_basis(nphot/2, modes/2);
	fock input_state(modes);
	for(int i = 0; i < modes/2; i++)
		input_state[i] = 1;
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
	using clock = chrono::high_resolution_clock;
	auto t1 = clock::now();
	auto out = C.output_state<execution::par>();
	auto t2 = clock::now();
	chrono::duration<double> sec = t2 - t1;
	//for(const auto &anc : ancilla_basis)
	//	cout << out.postselect(anc).normalize() << endl;
	cout << "output_state() execution time is " << sec.count() << " seconds with "
		 << omp_get_max_threads() << " threads." << endl;
	return EXIT_SUCCESS;
}

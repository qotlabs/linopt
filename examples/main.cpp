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
	Basis fullBasis = Basis(nphot/2, modes/2)*Basis(nphot/2, modes/2);
	Basis ancillaBasis(nphot/2, modes/2);
	Fock inputState(modes);
	for(int i = 0; i < modes/2; i++)
		inputState[i] = 1;
	Point a(modes*modes);
	random_device rand_dev;
	mt19937 randGen(rand_dev());
	uniform_real_distribution<Real> distribution(0., 1.);
	for(unsigned i = 0; i < a.size(); i++)
		a[i] = distribution(randGen);
	Circuit C;
	C.setUnitary(hurwitzParametrization(a));
	C.setInputState(inputState);
	C.setOutputBasis(fullBasis);
	using clock = chrono::high_resolution_clock;
	auto t1 = clock::now();
	auto out = C.outputState<execution::Par>();
	auto t2 = clock::now();
	chrono::duration<double> sec = t2 - t1;
	//for(const auto &anc : ancilla_basis)
	//	cout << out.postselect(anc).normalize() << endl;
	cout << "output_state() execution time is " << sec.count() << " seconds with "
		 << omp_get_max_threads() << " threads." << endl;
	return EXIT_SUCCESS;
}

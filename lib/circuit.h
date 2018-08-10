#ifndef CIRCUIT_H
#define CIRCUIT_H

#include "states.h"
#include "matrix.h"

namespace linopt
{

class circuit
{
private:
	matrix_type U;
	matrix_type Uin;
	bool uin_possibly_changed;
	matrix_type Uinout;
	fock _input_state;
	real_type input_prod_fact;
	basis _output_basis;
	matrix_type &prepare_uin(const matrix_type &u, const fock &fin);
	complex_type calc_fock_amp(const fock &fout);

public:
	circuit();
	matrix_type &unitary();
	const matrix_type &unitary() const;
	fock &input_state();
	const fock &input_state() const;
	basis &output_basis();
	const basis &output_basis() const;
	state output_state();
	void set_input_state(const fock &fin);
	void set_output_basis(const basis &bout);
	void set_unitary(const matrix_type &u);
};

}

#endif // CIRCUIT_H

#ifndef CIRCUIT_H
#define CIRCUIT_H

#include "states.h"
#include "matrix.h"

namespace linopt
{

class circuit
{
private:
    unitary_matrix U;
    unitary_matrix Uin;
    bool uin_possibly_changed;
    unitary_matrix Uinout;
    fock _input_state;
    real_type input_prod_fact;
    basis _output_basis;
    unitary_matrix &prepare_uin(const unitary_matrix &u, const fock &fin);
    complex_type calc_fock_amp(const fock &fout);

public:
	circuit();
    unitary_matrix &unitary();
    const unitary_matrix &unitary() const;
    fock &input_state();
    const fock &input_state() const;
    basis &output_basis();
    const basis &output_basis() const;
    state output_state();
    void set_input_state(const fock &fin);
    void set_output_basis(const basis &bout);
    void set_unitary(const unitary_matrix &u);
};

}

#endif // CIRCUIT_H

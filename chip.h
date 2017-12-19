#ifndef CHIP_H
#define CHIP_H

#include "states.h"
#include "matrix.h"

namespace linopt
{

class chip
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
    chip();
    unitary_matrix &unitary();
    const unitary_matrix &unitary() const;
    fock &input_state();
    const fock &input_state() const;
    basis &output_basis();
    const basis &output_basis() const;
    state output_state();
};

}

#endif // CHIP_H

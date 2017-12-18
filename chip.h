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
    fock input_state;
    real_type input_prod_fact;
    basis output_basis;
    unitary_matrix &prepare_uin(const unitary_matrix &u, const fock &fin);
    complex_type calc_fock_amp(const fock &fout);

public:
    chip();
    chip &set_unitary(const unitary_matrix &u);
    unitary_matrix &unitary();
    chip &set_input(const fock &f);
    fock &input();
    chip &set_basis(const basis &b);
    state output_state();
};

}

#endif // CHIP_H

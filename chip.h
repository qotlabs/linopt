#ifndef CHIP_H
#define CHIP_H

#include "states.h"
#include "matrix.h"

namespace linopt
{

class chip
{
private:
    unitary_matrix unitary;
    unitary_matrix uin_prepared;
    fock input_state;
    basis output_basis;
    unitary_matrix &prepare_uin(const unitary_matrix &U, const fock &fin);

public:
    chip();
    chip &set_unitary(const unitary_matrix &U);
    chip &set_input(const fock &f);
    chip &set_basis(const basis &b);
    state output_state() const;
};

}

#endif // CHIP_H

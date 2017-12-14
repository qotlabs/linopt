#include <iostream>
#include "chip.h"

using namespace linopt;

chip::chip():
    unitary(),
    uin_prepared(),
    input_state(),
    output_basis() {}

unitary_matrix &chip::prepare_uin(const unitary_matrix &U, const fock &fin)
{
    int tot = fin.total();
    int modes = fin.size();
    int k = 0;
    uin_prepared.resize(U.rows(), tot);
    for(int m = 0; m < modes; m++)
        for(int i = 0; i < fin[m]; i++)
            uin_prepared.col(k++) = U.col(m);
    std::cout << uin_prepared << std::endl;
    return uin_prepared;
}

chip &chip::set_unitary(const unitary_matrix &U)
{
    unitary = U;
    prepare_uin(unitary, input_state);
    return *this;
}

chip &chip::set_input(const fock &f)
{
    input_state = f;
    prepare_uin(unitary, input_state);
    return *this;
}

chip &chip::set_basis(const basis &b)
{
    output_basis = b;
    return *this;
}

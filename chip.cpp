#include "chip.h"
#include <functional>

using namespace linopt;

chip::chip():
    U(),
    Uin(),
    uin_possibly_changed(false),
    Uinout(),
    input_state(),
    output_basis() {}

unitary_matrix &chip::prepare_uin(const unitary_matrix &u, const fock &fin)
{
    int tot = fin.total();
    int modes = fin.size();
    int k = 0;
    Uin.resize(u.rows(), tot);
    for(int m = 0; m < modes; m++)
        for(int i = 0; i < fin[m]; i++)
            Uin.col(k++) = u.col(m);
    return Uin;
}

complex_type chip::calc_fock_amp(const fock &fout)
{
    int tot = fout.total();
    int modes = fout.size();
    Uinout.resize(tot, Uin.cols());
    int k = 0;
    for(int m = 0; m < modes; m++)
        for(int i = 0; i < fout[m]; i++)
            Uinout.row(k++) = Uin.row(m);
    complex_type perm = permanent(Uinout);
    perm /= std::sqrt(fout.prod_fact() * input_state.prod_fact());
    return perm;
}

chip &chip::set_unitary(const unitary_matrix &u)
{
    U = u;
    uin_possibly_changed = true;
    return *this;
}

unitary_matrix &chip::unitary()
{
    uin_possibly_changed = true;
    return U;
}

chip &chip::set_input(const fock &f)
{
    input_state = f;
    uin_possibly_changed = true;
    return *this;
}

fock &chip::input()
{
    uin_possibly_changed = true;
    return input_state;
}

chip &chip::set_basis(const basis &b)
{
    output_basis = b;
    return *this;
}

state chip::output_state()
{
    if(uin_possibly_changed)
    {
        prepare_uin(U, input_state);
        uin_possibly_changed = false;
    }
    return output_basis.apply_func(std::bind(&chip::calc_fock_amp,
                                             this, std::placeholders::_1));
}

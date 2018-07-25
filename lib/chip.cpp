#include "chip.h"
#include <functional>

using namespace linopt;

chip::chip():
    uin_possibly_changed(false) {}

unitary_matrix &chip::prepare_uin(const unitary_matrix &u, const fock &fin)
{
    int tot = fin.total();
    int modes = fin.size();
    int k = 0;
    Uin.resize(u.rows(), tot);
    for(int m = 0; m < modes; m++)
        for(int i = 0; i < fin[m]; i++)
            Uin.col(k++) = u.col(m);
    input_prod_fact = fin.prod_fact();
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
    perm /= std::sqrt(fout.prod_fact() * input_prod_fact);
    return perm;
}

unitary_matrix &chip::unitary()
{
    uin_possibly_changed = true;
    return U;
}

const unitary_matrix &chip::unitary() const
{
    return U;
}

fock &chip::input_state()
{
    uin_possibly_changed = true;
    return _input_state;
}

const fock &chip::input_state() const
{
    return _input_state;
}

basis &chip::output_basis()
{
    return _output_basis;
}

const basis &chip::output_basis() const
{
    return _output_basis;
}

state chip::output_state()
{
    if(uin_possibly_changed)
    {
        prepare_uin(U, _input_state);
        uin_possibly_changed = false;
    }
    return _output_basis.apply_func(std::bind(&chip::calc_fock_amp,
                                             this, std::placeholders::_1));
}

void chip::set_input_state(const fock &fin)
{
    _input_state = fin;
}

void chip::set_output_basis(const basis &bout)
{
    _output_basis = bout;
}

void chip::set_unitary(const unitary_matrix &u)
{
    U = u;
}
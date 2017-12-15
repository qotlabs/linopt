#ifndef COST_FUNCTOR_H
#define COST_FUNCTOR_H

#include "chip.h"
#include <vector>

namespace linopt
{

class cost_functor
{
protected:
    std::vector<state> target_states;
    basis ancilla_basis;
    chip C;
public:
    cost_functor(const basis &full_basis,
                 const basis &ancilla_basis,
                 const fock &input_state,
                 const std::vector<state> &target_states);
    real_type operator()(const unitary_matrix::angles &a);
};

class stanisic_functor: public cost_functor
{
public:
    using cost_functor::cost_functor;
    real_type operator()(const unitary_matrix::angles &a);
};

}

#endif // COST_FUNCTOR_H

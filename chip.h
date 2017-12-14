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
public:
    chip();
    chip &setUnitary(const unitary_matrix &U);

};

}

#endif // CHIP_H

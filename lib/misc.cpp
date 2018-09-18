#include "misc.h"

using namespace linopt;

std::ostream& linopt::print_complex(std::ostream &stream, const complex_type &x)
{
	const real_type re = x.real();
	const real_type im = x.imag();
	if(re < 0)
		stream << "-";
	stream << re;
	if(im < 0)
		stream << "-";
	else
		stream << "+";
	stream << im << "j";
	return stream;
}

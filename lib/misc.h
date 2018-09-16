#ifndef MISC_H
#define MISC_H

#include <cmath>
#include <ostream>

namespace linopt
{

template<typename T>
inline T mod(T x, T a)
{
	x = std::fmod(x, a);
	return (x < 0) ? x + a : x;
}

template<typename T>
std::ostream& print_array(std::ostream &stream, const T &a,
						  const char *b1 = "{",
						  const char *delim = ", ",
						  const char *b2 = "}")
{
	if(a.empty())
	{
		stream << b1 << b2;
	}
	else
	{
		stream << b1;
		for(auto iter = a.begin(); iter != --a.end(); iter++)
			stream << *iter << delim;
		stream << *(--a.end()) << b2;
	}
	return stream;
}

}

#endif // MISC_H

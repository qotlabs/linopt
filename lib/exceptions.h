#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <string>
#include <stdexcept>

namespace linopt
{

inline std::string format_error_message(std::string msg, const char *func)
{
	return std::string() + "linopt::" + func + "(): " + msg;
}

#define ERROR_MSG(...) format_error_message(__VA_ARGS__, __func__)

class general_error: public std::logic_error
{
	using logic_error::logic_error;
};

class wrong_size: public general_error
{
	using general_error::general_error;
};

}

#endif // EXCEPTIONS_H

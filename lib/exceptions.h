/* Copyright © 2018, Quantum Optical Technologies Laboratories
 * <https://www.qotlabs.org/en/>
 * Contributed by: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
 *                 Dyakonov Ivan <iv.dyakonov@physics.msu.ru>
 *
 * This file is part of Linopt.
 *
 * Linopt is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Linopt is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Linopt. If not, see <https://www.gnu.org/licenses/>.
 */

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

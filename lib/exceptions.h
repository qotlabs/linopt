/* Copyright © 2018-2020, Quantum Optical Technologies Laboratories
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

/** @defgroup exceptions Exceptions
  * @brief Exception types.
  */

#ifndef _LINOPT_EXCEPTIONS_H
#define _LINOPT_EXCEPTIONS_H

#include <string>
#include <stdexcept>

namespace linopt
{

inline std::string formatErrorMessage(std::string msg, const char *func)
{
	return std::string() + "linopt::" + func + "(): " + msg;
}

#define ERROR_MSG(...) formatErrorMessage(__VA_ARGS__, __func__)

/** @ingroup exceptions
 * @brief The base exception class.
 *
 * All other exceptions inherits from it. Therefore, you can catch objects of
 * this type to handle all Linopt exceptions.
 */
class GeneralError: public std::logic_error
{
	using logic_error::logic_error;
};

/** @ingroup exceptions
 * @brief The object of this type is thrown when an object of improper size
 * is encountered.
 */
class WrongSize: public GeneralError
{
	using GeneralError::GeneralError;
};

/** @ingroup exceptions
 * @brief The object of this type is thrown when a unitary matrix is expected,
 * but it is not.
 */
class NotUnitary : public GeneralError
{
	using GeneralError::GeneralError;
};

class NotSupported: public GeneralError
{
	using GeneralError::GeneralError;
};

} // Namespace linopt

#endif // _LINOPT_EXCEPTIONS_H

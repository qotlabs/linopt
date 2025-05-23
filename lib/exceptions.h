// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

/** @defgroup exceptions Exceptions
  * @brief Exception types.
  */

#pragma once

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

} // namespace linopt

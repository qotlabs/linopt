// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

#pragma once

#include <complex>

namespace linopt
{

using Real = double;
using Complex = std::complex<Real>;

/** @ingroup matrix
 * @brief Default precision for numeric comparison operations.
 */
constexpr Real defaultEpsilon = 1e-15;

} // namespace linopt

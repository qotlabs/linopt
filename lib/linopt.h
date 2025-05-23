// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

/** @mainpage
 * Linopt is a numerical library designed for linear-optical quantum system
 * simulation. The core functions are implemented using C++ to enable efficient
 * performance. The classes  and functions are exposed to Python which enables
 * simple and easy use of the library.
 *
 * The package is designed to simulate the evolution of bosonic states encoded
 * in occupation number basis through a discrete unitary network. The generic
 * scheme of the system is depicted below.
 *
 * @image html lop_gate_general.png "Linear optical gate" width=400em
 *
 * The package implements methods for computing probabilities of finding the
 * output state in the fock basis given the defined input state.
 * The implementation allows easy manipulation with the bosonic fock states,
 * the basis state sets and the unitary networks. The package is primarily
 * developed to enable efficient design of linear-optical interferometer for
 * quantum computing purposes. It enables the user to easily compute the
 * postselected output state of a given linear-optical interferometer or deliver
 * the configuration of the interferometer given a precomputed unitary matrix.
 */

#pragma once

#include "types.h"
#include "matrix.h"
#include "states.h"
#include "circuit.h"
#include "circuit_design.h"

/**
 * @brief The main namespace containing all library classes, functions, etc.
 */
namespace linopt
{

/**
  * @brief Library version specification.
  */
struct Version
{
	static constexpr auto string = "0.4.0";
	static constexpr int major = 0;
	static constexpr int minor = 4;
	static constexpr int patch = 0;
};

} // namespace linopt

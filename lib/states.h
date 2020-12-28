/* Copyright © 2018, 2019, Quantum Optical Technologies Laboratories
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

#ifndef STATES_H
#define STATES_H

#include <vector>
#include <set>
#include <map>
#include <ostream>
#include <initializer_list>
#include <functional>
#include "types.h"

namespace linopt
{

class Fock;
class Basis;
class State;

/** @ingroup states
 * @brief A typedef of a function taking a Fock state as an argument and
 * returning a corresponding `Complex` amplitude.
 */
using FockAmpFunction = std::function<Complex(const Fock&)>;

/** @ingroup states
 * @brief The class representing a Fock state.
 *
 * A Fock state is a state with exact number of photons in each mode. A Fock
 * state is represented as a list of nonnegative integers:
 * @f[ | n_0, \dots, n_m \rangle, @f]
 * where @f$ n_i @f$ is called an _occupation number_ of an @f$ i @f$-th mode.
 *
 * The interface of this class is very similar to `std::vector`.
 */
class Fock: private std::vector<int>
{
	using Base = std::vector<int>;

public:
	/// Convenience typedef to std::vector.
	using Vector = std::vector<int>; // In principle this can be different from Base
	using Value = int;
	using value_type = Value;
	using Reference = Value&;
	using reference = Reference;

	using Base::Base;
	using Base::operator=;

	using Iterator = Base::iterator;
	using ConstIterator = Base::const_iterator;
	using ReverseIterator = Base::reverse_iterator;
	using ConstReverseIterator = Base::const_reverse_iterator;

	using Base::begin;
	using Base::end;
	using Base::rbegin;
	using Base::rend;
	using Base::cbegin;
	using Base::cend;
	using Base::crbegin;
	using Base::crend;

	using Base::front;
	using Base::back;
	using Base::operator[];

	/// Checks whether a Fock is empty, i.e. contains zero modes.
	using Base::empty;
	/// Returns the number of modes in a Fock state.
	int size() const { return static_cast<int>(Base::size()); }
	using Base::resize;
	using Base::assign;
	void pushBack(Value n) { Base::push_back(n); }
	void popBack() { Base::pop_back(); }
	using Base::insert;
	using Base::erase;
	/// Clears a Fock state.
	using Base::clear;

	/// Default constructor.
	Fock() = default;
	/// Construct a Fock state from a `fock::vector_class`.
	Fock(const Vector &v): Base(v) {}
	int total() const;
	Real prodFact() const;
	Fock operator*(const Fock &f) const;
	Fock &operator*=(const Fock &f);
	Fock operator+(const Fock &f) const;
	Fock &operator+=(const Fock &f);
};

/** @ingroup states
 * @brief Tests whether two Fock states are equal.
 */
inline bool operator==(const Fock &f, const Fock &g)
{
	return f.size() == g.size() &&
			std::equal(f.begin(), f.end(), g.begin());
}

/** @ingroup states
 * @brief Compares two Fock states in lexicographic order.
 *
 * @return
 * `true` if `f` is less than `g`.
 */
inline bool operator<(const Fock &f, const Fock &g)
{
	return std::lexicographical_compare(f.begin(), f.end(),
										g.begin(), g.end());
}

/** @ingroup states
 * @brief Tests whether two Fock states differ.
 */
inline bool operator!=(const Fock &f, const Fock &g)
{
	return !(f == g);
}

/** @ingroup states
 * @brief Compares two Fock states in lexicographic order.
 *
 * @return
 * `true` if `f` is greater than `g`.
 */
inline bool operator>(const Fock &f, const Fock &g)
{
	return g < f;
}

/** @ingroup states
 * @brief Compares two Fock states in lexicographic order.
 *
 * @return
 * `true` if `f` is less or equal than `g`.
 */
inline bool operator<=(const Fock &f, const Fock &g)
{
	return !(f > g);
}

/** @ingroup states
 * @brief Compares two Fock states in lexicographic order.
 *
 * @return
 * `true` if `f` is greater or equal than `g`.
 */
inline bool operator>=(const Fock &f, const Fock &g)
{
	return !(f < g);
}

/** @ingroup states
 * @brief The class representing a collection of Fock states.
 *
 * The interface of this class is very similar to `std::set<fock>`.
 */
class Basis: private std::set<Fock>
{
	using Base = std::set<Fock>;

public:
	/// Convenience typedef to std::set.
	using Set = std::set<Fock>; // In principle this can be different from Base
	using Value = Fock;
	using value_type = Value;
	using Reference = Value&;
	using reference = Reference;

	using Iterator = Base::iterator;
	using ConstIterator = Base::const_iterator;
	using ReverseIterator = Base::reverse_iterator;
	using ConstReverseIterator = Base::const_reverse_iterator;

	using Base::begin;
	using Base::end;
	using Base::rbegin;
	using Base::rend;
	using Base::cbegin;
	using Base::cend;
	using Base::crbegin;
	using Base::crend;

	using Base::empty;
	/// Returns the number of Fock states in the basis.
	int size() const { return static_cast<int>(Base::size()); }
	using Base::insert;
	using Base::erase;
	using Base::clear;
	using Base::find;

	using Base::operator=;
	Basis(): Base() {}
	Basis(const Set &s): Base(s) {}
	Basis(std::initializer_list<Fock> il): Base(il) {}
	explicit Basis(int nphot, int modes);

	Basis operator+(const Basis &b) const;
	Basis &operator+=(const Basis &b);
	Basis operator*(const Basis &b) const;
	Basis &operator*=(const Basis &b);
	Basis &generateBasis(const int nphot, const int modes, const Fock &head = Fock());
	Basis postselect(const Fock &ancilla) const;
	template<typename ExecPolicy = execution::Seq>
	State applyFunction(const FockAmpFunction &f) const;
};

/** @ingroup states
 * @brief The class representing a linear optical state.
 *
 * The interface of this class is very similar to `std::map<Fock, Complex>`.
 */
class State: private std::map<Fock, Complex>
{
	using Base = std::map<Fock, Complex>;

public:
	using Element = std::pair<Fock, Complex>;
	using Map = std::map<Fock, Complex>; // In principle this can be different from Base
	using Value = Complex;
	using value_type = Value;
	using Reference = Value&;
	using reference = Reference;

	using Iterator = Base::iterator;
	using ConstIterator = Base::const_iterator;
	using ReverseIterator = Base::reverse_iterator;
	using ConstReverseIterator = Base::const_reverse_iterator;

	using Base::begin;
	using Base::end;
	using Base::rbegin;
	using Base::rend;
	using Base::cbegin;
	using Base::cend;
	using Base::crbegin;
	using Base::crend;

	using Base::empty;
	int size() const { return static_cast<int>(Base::size()); }
	using Base::operator[];
	using Base::insert;
	using Base::erase;
	using Base::clear;
	using Base::find;

	using Base::operator=;
	/// Default constructor.
	State(): Base() {}
	State(const Map &m): Base(m) {}
	/// Constructs a `state` from `fock` with unit amplitude.
	State(const Fock &f): Base() { (*this)[f] = 1; }
	State(const Basis &b): Base() { setBasis(b); }

	State operator+(const State &s) const;
	State &operator+=(const State &s);
	/// Subtracts two states.
	State operator-(const State &s) const { return *this + (-s); }
	State &operator-=(const State &s);
	State operator*(const State &s) const;
	State &operator*=(const State &s);
	State operator-() const;
	State operator*(Complex x) const;
	State &operator*=(Complex x);
	State operator/(Complex x) const;
	State &operator/=(Complex x);
	Real norm() const;
	State &normalize();
	Complex dot(const State &s) const;
	State postselect(const Fock &ancilla) const;
	std::map<Fock, State> postselect(int modes) const;
	std::map<Fock, State> postselect(const Basis &b) const;
	Basis getBasis() const;
	void setBasis(const Basis &b);
	std::vector<State::Value> getAmplitudes() const;
	template<typename ExecPolicy = execution::Seq>
	void setAmplitudes(const std::vector<Complex> &amps);
	template<typename ExecPolicy = execution::Seq>
	void setAmplitudes(const FockAmpFunction &f);
};

State operator*(Complex x, const State &s);
Complex dot(const State &a, const State &b);

} // Namespace linopt

std::ostream &operator<<(std::ostream &stream, const linopt::Fock &f);
std::ostream &operator<<(std::ostream &stream, const linopt::Basis &b);
std::ostream &operator<<(std::ostream &stream, const linopt::State::Element &e);
std::ostream &operator<<(std::ostream &stream, const linopt::State &s);

#endif // STATES_H

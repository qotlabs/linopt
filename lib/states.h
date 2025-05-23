// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright (c) 2018-2025, Quantum Optical Technologies Laboratories
// SPDX-FileContributor: Struchalin Gleb <struchalin.gleb@physics.msu.ru>
// SPDX-FileContributor: Dyakonov Ivan <iv.dyakonov@physics.msu.ru>

#pragma once

#include "types.h"
#include <vector>
#include <set>
#include <map>
#include <ostream>
#include <initializer_list>
#include <functional>
#include <execution>

/** @defgroup states States
 * @brief Linear optics states.
 */

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

	/// Check whether a Fock is empty, i.e. contains zero modes.
	using Base::empty;

	/// Return the number of modes in a Fock state.
	int size() const noexcept { return static_cast<int>(Base::size()); }

	using Base::resize;
	using Base::assign;
	void pushBack(Value n) { Base::push_back(n); }
	void popBack() { Base::pop_back(); }
	using Base::insert;
	using Base::erase;

	/// Clear a Fock state.
	using Base::clear;

	/// Default constructor.
	Fock() = default;

	/// Construct a Fock state from a `fock::vector_class`.
	Fock(const Vector &v): Base(v) {}

	/**
	 * @brief Return the total number of photons in all modes.
	 */
	int total() const;

	/**
	 * @brief Return a product of factorials of occupation numbers.
	 */
	Real prodFact() const;

	/**
	 * @brief Return a tensor product of `*this` and `f`.
	 */
	Fock operator*(const Fock &f) const;

	/**
	 * @brief Effectively equivalent to `*this = (*this) * f`.
	 */
	Fock &operator*=(const Fock &f);

	/**
	 * @brief Return a sum of two Fock states (elementwise addition of
	 * corresponding occupation numbers).
	 *
	 * @throw
	 * If `*this` and `f` have different sizes then `WrongSize` is thrown.
	 */
	Fock operator+(const Fock &f) const;

	/**
	 * @brief Effectively equivalent to `*this = (*this) + f`.
	 */
	Fock &operator+=(const Fock &f);
};

/** @ingroup states
 * @brief Check whether two Fock states are equal.
 */
inline bool operator==(const Fock &f, const Fock &g)
{
	return f.size() == g.size() &&
			std::equal(f.begin(), f.end(), g.begin());
}

/** @ingroup states
 * @brief Compare two Fock states in lexicographic order.
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
 * @brief Test whether two Fock states differ.
 */
inline bool operator!=(const Fock &f, const Fock &g)
{
	return !(f == g);
}

/** @ingroup states
 * @brief Compare two Fock states in lexicographic order.
 *
 * @return
 * `true` if `f` is greater than `g`.
 */
inline bool operator>(const Fock &f, const Fock &g)
{
	return g < f;
}

/** @ingroup states
 * @brief Compare two Fock states in lexicographic order.
 *
 * @return
 * `true` if `f` is less or equal than `g`.
 */
inline bool operator<=(const Fock &f, const Fock &g)
{
	return !(f > g);
}

/** @ingroup states
 * @brief Compare two Fock states in lexicographic order.
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

	/// Return the number of Fock states in the basis.
	int size() const noexcept { return static_cast<int>(Base::size()); }

	using Base::insert;
	using Base::erase;
	using Base::clear;
	using Base::find;

	/// Check whether two bases are equal.
	bool operator ==(const Basis &b) const
	{
		return static_cast<const Base&>(*this) == static_cast<const Base&>(b);
	}

	/// Check whether two bases differ.
	bool operator !=(const Basis &b) const
	{
		return static_cast<const Base&>(*this) != static_cast<const Base&>(b);
	}

	using Base::operator=;
	Basis(): Base() {}
	Basis(const Set &s): Base(s) {}
	Basis(std::initializer_list<Fock> il): Base(il) {}

	/**
	 * @brief Construct a basis of all possible Fock states with `modes` modes
	 * and containing `nphot` photons.
	 */
	explicit Basis(int nphot, int modes);

	/**
	 * @brief Return a basis which is a union of Fock states from both `*this`
	 * and `b`.
	 */
	Basis operator+(const Basis &b) const;

	/**
	 * @brief Effectively equivalent to `*this = (*this) + b`.
	 */
	Basis &operator+=(const Basis &b);

	/**
	 * @brief Calculate a tensor product of two bases.
	 *
	 * Return a basis consisting of all possible elementwise tensor products
	 * of elements of `*this` and `b`.
	 */
	Basis operator*(const Basis &b) const;

	/**
	 * @brief Effectively equivalent to `*this = (*this) * b`.
	 */
	Basis &operator*=(const Basis &b);

	/**
	 * @brief Generate a basis of all possible Fock states with `modes` modes
	 * and containing `nphot` photons.
	 *
	 * @param[in] nphot -- number of photons in each Fock state.
	 * @param[in] modes -- number of modes in each Fock state.
	 * @deprecated
	 * @param[in] head -- intended for internal usage. Normally empty Fock state
	 * should be passed.
	 *
	 * @note
	 * This funtion only appends elements to `*this` and never removes them.
	 * Therefore, if you want to freshly generate a basis, you should call
	 * `clear()` first.
	 */
	Basis &generateBasis(int nphot, int modes, const Fock &head = Fock());

	/**
	 * @brief Return a postselected basis after observing ancilla.
	 *
	 * @param[in] ancilla -- ancilla's Fock state for postselection.
	 * @return A basis after postselection.
	 *
	 * Postselection of a basis @f$ B @f$ with the ancilla state
	 * @f$ | \mathrm{anc} \rangle = | a_0, a_2, \dots, a_A \rangle @f$
	 * picks only Fock states from @f$ B @f$ with the heading (first @f$ A @f$
	 * occupation numbers) which equals to @f$ | \mathrm{anc} \rangle @f$.
	 * The headings do not get to the new constructed basis.
	 * Note that currently there is no possibility to specify a position of
	 * ancilla modes.
	 *
	 * For example, postselection of a basis
	 * @f[ \{ | 0000001 \rangle, | 1234567 \rangle, | 1239999 \rangle \} @f]
	 * with the ancilla @f$ | 123 \rangle @f$ results in the basis
	 * @f[ \{ | 4567 \rangle, | 9999 \rangle \}. @f]
	 *
	 */
	Basis postselect(const Fock &ancilla) const;

	/**
	 * @brief Construct a state from the basis using function `f`.
	 *
	 * @param[in] f -- function to apply.
	 * @return Constructed state.
	 *
	 * Apply a function `f` to all Fock states of `*this` to compute a
	 * corresponding amplitude of a resulting state.
	 */
	template<typename ExecPolicy = std::execution::sequenced_policy>
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

	/// Construct this state from `Fock` with unit amplitude.
	State(const Fock &f): Base() { (*this)[f] = 1; }

	State(const Basis &b): Base() { setBasis(b); }

	/// Check whether two states are equal.
	bool operator ==(const State &s) const
	{
		return static_cast<const Base&>(*this) == static_cast<const Base&>(s);
	}

	/// Check whether two states differ.
	bool operator !=(const State &s) const
	{
		return static_cast<const Base&>(*this) != static_cast<const Base&>(s);
	}

	/**
	 * @brief Add two states, i.e., calculate their superposition.
	 */
	State operator+(const State &s) const;

	/**
	 * @brief Effectively equivalent to `*this = (*this) + s`.
	 */
	State &operator+=(const State &s);

	/// Subtract two states.
	State operator-(const State &s) const { return *this + (-s); }

	/**
	 * @brief Effectively equivalent to `*this = (*this) - s`.
	 */
	State &operator-=(const State &s);

	/**
	 * @brief Return a tensor product of two states.
	 */
	State operator*(const State &s) const;

	/**
	 * @brief Effectively equivalent to `*this = (*this) * s`.
	 */
	State &operator*=(const State &s);

	/**
	 * @brief Negate amplitudes of the state.
	 */
	State operator-() const;

	/**
	 * @brief Multiply a state by a complex number.
	 */
	State operator*(Complex x) const;

	/**
	 * @brief Effectively equivalent to `*this = (*this) * x`.
	 */
	State &operator*=(Complex x);

	/**
	 * @brief Divide a state by a complex number.
	 */
	State operator/(Complex x) const;

	/**
	 * @brief Effectively equivalent to `*this = (*this) / x`.
	 */
	State &operator/=(Complex x);

	/**
	 * @brief Return norm of the state.
	 */
	Real norm() const;

	/**
	 * @brief Normalize the state to have unit norm.
	 */
	State &normalize();

	/**
	 * @brief Calculates a dot (scalar) product.
	 */
	Complex dot(const State &s) const;

	State postselect(const Fock &ancilla) const;

	std::map<Fock, State> postselect(int modes) const;

	std::map<Fock, State> postselect(const Basis &b) const;

	Basis getBasis() const;

	void setBasis(const Basis &b);

	std::vector<State::Value> getAmplitudes() const;

	template<typename ExecPolicy = std::execution::sequenced_policy>
	void setAmplitudes(const std::vector<Complex> &amps);

	template<typename ExecPolicy = std::execution::sequenced_policy>
	void setAmplitudes(const FockAmpFunction &f);
};

/** @ingroup states
 * @brief Multiply a state by a complex number.
 */
State operator*(Complex x, const State &s);

/** @ingroup states
 * @brief Calculate a dot (scalar) product.
 */
static inline Complex dot(const State &a, const State &b) { return a.dot(b); }

} // Namespace linopt

std::ostream &operator<<(std::ostream &stream, const linopt::Fock &f);
std::ostream &operator<<(std::ostream &stream, const linopt::Basis &b);
std::ostream &operator<<(std::ostream &stream, const linopt::State::Element &e);
std::ostream &operator<<(std::ostream &stream, const linopt::State &s);

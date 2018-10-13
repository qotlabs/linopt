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

#ifndef STATES_H
#define STATES_H

#include <vector>
#include <set>
#include <map>
#include <ostream>
#include <initializer_list>
#include "types.h"

namespace linopt
{

class fock : private std::vector<int>
{
	typedef std::vector<int> base_class;
public:
	typedef std::vector<int> vector_class; //In principle this can be different from base_class
	typedef int value_type;
	typedef int& reference;

	using base_class::base_class;
	using base_class::operator=;

	using base_class::iterator;
	using base_class::const_iterator;
	using base_class::reverse_iterator;
	using base_class::const_reverse_iterator;

	using base_class::begin;
	using base_class::end;
	using base_class::rbegin;
	using base_class::rend;
	using base_class::cbegin;
	using base_class::cend;
	using base_class::crbegin;
	using base_class::crend;

	using base_class::front;
	using base_class::back;
	using base_class::operator[];

	using base_class::empty;
	int size() const { return base_class::size(); }
	using base_class::resize;
	using base_class::assign;
	using base_class::push_back;
	using base_class::pop_back;
	using base_class::insert;
	using base_class::erase;
	using base_class::clear;

	fock(const vector_class &v): base_class(v) {}
	int total() const;
	real_type prod_fact() const;
	fock operator*(const fock &f) const;
	fock &operator*=(const fock &f);
	fock operator+(const fock &f) const;
	fock &operator+=(const fock &f);
};

inline bool operator==(const fock &f, const fock &g)
{
	return f.size() == g.size() &&
			std::equal(f.begin(), f.end(), g.begin());
}

inline bool operator<(const fock &f, const fock &g)
{
	return std::lexicographical_compare(f.begin(), f.end(),
										g.begin(), g.end());
}

inline bool operator!=(const fock &f, const fock &g)
{
	return !(f == g);
}

inline bool operator>(const fock &f, const fock &g)
{
	return g < f;
}

inline bool operator<=(const fock &f, const fock &g)
{
	return !(f > g);
}

inline bool operator>=(const fock &f, const fock &g)
{
	return !(f < g);
}

class basis : private std::set<fock>
{
	typedef std::set<fock> base_class;
public:
	typedef std::set<fock> set_class; //In principle this can be different from base_class
	typedef fock value_type;
	typedef fock& reference;

	using base_class::iterator;
	using base_class::const_iterator;
	using base_class::reverse_iterator;
	using base_class::const_reverse_iterator;

	using base_class::begin;
	using base_class::end;
	using base_class::rbegin;
	using base_class::rend;
	using base_class::cbegin;
	using base_class::cend;
	using base_class::crbegin;
	using base_class::crend;

	using base_class::empty;
	int size() const { return base_class::size(); }
	using base_class::insert;
	using base_class::erase;
	using base_class::clear;
	using base_class::find;

	using base_class::operator=;
	basis(): base_class() {}
	basis(const set_class &s): base_class(s) {}
	basis(std::initializer_list<fock> il): base_class(il) {}
	explicit basis(int nphot, int modes);

	basis operator+(const basis &b) const;
	basis &operator+=(const basis &b);
	basis operator*(const basis &b) const;
	basis &operator*=(const basis &b);
	basis &generate_basis(const int nphot, const int modes, const fock &head = fock());
	basis postselect(const fock &ancilla) const;
	state apply_func(const basis_func &f) const;
};

class state : private std::map<fock, complex_type>
{
	typedef std::map<fock, complex_type> base_class;
public:
	typedef std::pair<fock, complex_type> element;
	typedef std::map<fock, complex_type> map_class; //In principle this can be different from base_class
	typedef complex_type value_type;
	typedef complex_type& reference;

	using base_class::iterator;
	using base_class::const_iterator;
	using base_class::reverse_iterator;
	using base_class::const_reverse_iterator;

	using base_class::begin;
	using base_class::end;
	using base_class::rbegin;
	using base_class::rend;
	using base_class::cbegin;
	using base_class::cend;
	using base_class::crbegin;
	using base_class::crend;

	using base_class::empty;
	int size() const { return base_class::size(); }
	using base_class::operator[];
	using base_class::insert;
	using base_class::erase;
	using base_class::clear;
	using base_class::find;

	using base_class::operator=;
	state(): base_class() {}
	state(const map_class &m): base_class(m) {}
	state(const fock &f): base_class() { (*this)[f] = 1; }

	state operator+(const state &s) const;
	state &operator+=(const state &s);
	state operator-(const state &s) const { return *this + (-s); }
	state &operator-=(const state &s);
	state operator*(const state &s) const;
	state &operator*=(const state &s);
	state operator-() const;
	state operator*(complex_type x) const;
	state &operator*=(complex_type x);
	state operator/(complex_type x) const;
	state &operator/=(complex_type x);
	real_type norm() const;
	state &normalize();
	complex_type dot(const state &s) const;
	state postselect(const fock &ancilla) const;
	std::map<fock, state> postselect(int modes) const;
	std::map<fock, state> postselect(const basis &b) const;
	basis get_basis() const;
	std::vector<state::value_type> get_amplitudes() const;
};

state operator*(complex_type x, const state &s);
complex_type dot(const state &a, const state &b);

}

std::ostream &operator<<(std::ostream &stream, const linopt::fock &f);
std::ostream &operator<<(std::ostream &stream, const linopt::basis &b);
std::ostream &operator<<(std::ostream &stream, const linopt::state::element &e);
std::ostream &operator<<(std::ostream &stream, const linopt::state &s);

#endif // STATES_H

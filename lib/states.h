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
public:
	typedef std::vector<int> vector;

	using vector::vector;
	using vector::operator=;

	using vector::iterator;
	using vector::const_iterator;
	using vector::reverse_iterator;
	using vector::const_reverse_iterator;
	using vector::pointer;
	using vector::reference;

	using vector::begin;
	using vector::end;
	using vector::rbegin;
	using vector::rend;
	using vector::cbegin;
	using vector::cend;
	using vector::crbegin;
	using vector::crend;

	using vector::size;
	using vector::resize;
	using vector::empty;
	using vector::operator[];
	using vector::assign;
	using vector::push_back;
	using vector::pop_back;
	using vector::insert;
	using vector::erase;
	using vector::clear;

	fock(const vector &v): vector(v) {}
	int total() const;
	real_type prod_fact() const;
	fock operator*(const fock &f) const;
	fock &operator*=(const fock &f);
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

class basis : public std::set<fock>
{
	typedef std::set<fock> base_class;
public:
	basis(): base_class() {}
	basis(const basis &b): base_class(b) {}
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

class state : public std::map<fock, complex_type>
{
public:
	state();
	state(const state &s);
	state operator+(const state &s) const;
	state &operator+=(const state &s);
	state operator-(const state &s) const {return *this + (-s);}
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
	basis get_basis() const;
};

state operator*(complex_type x, const state &s);
complex_type dot(const state &a, const state &b);

}

std::ostream &operator<<(std::ostream &stream, const linopt::fock &f);
std::ostream &operator<<(std::ostream &stream, const linopt::basis &b);
std::ostream &operator<<(std::ostream &stream, const linopt::state_element &e);
std::ostream &operator<<(std::ostream &stream, const linopt::state &s);

#endif // STATES_H

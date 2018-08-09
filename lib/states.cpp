#include "states.h"

#include <algorithm>

using namespace linopt;

int fock::total() const
{
	int tot = 0;
	for(auto iter = begin(); iter != end(); iter++)
		tot += *iter;
	return tot;
}

static const int factorial_precomputed[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880,\
							   3628800, 39916800, 479001600};

static inline real_type factorial(int n)
{
	return (n <= 12) ? factorial_precomputed[n] : std::tgamma(n+1);
}

real_type fock::prod_fact() const
{
	real_type p = 1.;
	for(auto iter = begin(); iter != end(); iter++)
		p *= factorial(*iter);
	return p;
}

fock fock::operator*(const fock &f) const
{
	fock newf = *this;
	return newf *= f;
}

fock &fock::operator*=(const fock &f)
{
	insert(end(), f.begin(), f.end());
	return *this;
}

basis::basis():
	std::set<fock>() {}

basis::basis(const basis &b):
	std::set<fock>(b) {}

basis::basis(std::initializer_list<fock> il):
	std::set<fock>(il) {}

basis::basis(int nphot, int modes, const fock &head):
	basis()
{
	generate_basis(nphot, modes, head);
}

basis basis::operator+(const basis &b) const
{
	basis newb = *this;
	return newb += b;
}

basis &basis::operator+=(const basis &b)
{
	insert(b.begin(), b.end());
	return *this;
}

basis basis::operator*(const basis &b) const
{
	basis newb;
	for(auto f1 = begin(); f1 != end(); f1++)
		for(auto f2 = b.begin(); f2 != b.end(); f2++)
			newb.insert(newb.end(), (*f1) * (*f2));
	return newb;
}

basis &basis::operator*=(const basis &b)
{
	*this = (*this) * b;
	return *this;
}

basis &basis::generate_basis(const int nphot, const int modes, const fock &head)
{
	if(static_cast<int>(head.size()) == modes)
	{
		if(nphot == 0)
			this->insert(this->end(), head);
		return *this;
	}
	for(int i = 0; i <= nphot; i++)
	{
		fock f(head);
		f.push_back(i);
		generate_basis(nphot - i, modes, f);
	}
	return *this;
}

basis basis::postselect(const fock &ancilla) const
{
	basis b;
	for(auto fp = begin(); fp != end(); fp++)
	{
		if(std::equal(ancilla.rbegin(), ancilla.rend(), fp->rbegin()))
			b.insert(b.end(),
					 fock(fp->begin(), fp->end() - ancilla.size()));
	}
	return b;
}

state basis::apply_func(const basis_func &f) const
{
	state s;
	for(auto iter = begin(); iter != end(); iter++)
		s.insert(s.end(), state_element(*iter, f(*iter)));
	return s;
}

state::state():
	std::map<fock, complex_type>() {}

state::state(const state &s):
	std::map<fock, complex_type>(s) {}

state state::operator+(const state &s) const
{
	state snew;
	auto a = this->begin();
	auto b = s.begin();
	while(a != this->end() && b != s.end())
	{
		if(a->first < b->first)
		{
			snew.insert(snew.end(), *a);
			a++;
		}
		else if(a->first == b->first)
		{
			snew.insert(snew.end(), state_element(a->first, a->second + b->second));
			a++;
			b++;
		}
		else
		{
			snew.insert(snew.end(), *b);
			b++;
		}
	}
	if(a == this->end())
		snew.insert(b, s.end());
	else
		snew.insert(a, this->end());
	return snew;
}

state &state::operator+=(const state &s)
{
	for(auto iter = s.begin(); iter != s.end(); iter++)
		(*this)[iter->first] += iter->second;
	return *this;
}

state &state::operator-=(const state &s)
{
	for(auto iter = s.begin(); iter != s.end(); iter++)
		(*this)[iter->first] -= iter->second;
	return *this;
}

state state::operator*(const state &s) const
{
	state snew;
	for(auto a = begin(); a != end(); a++)
		for(auto b = s.begin(); b != s.end(); b++)
			snew.insert(snew.end(), state_element(a->first * b->first, a->second * b->second));
	return snew;
}

state &state::operator*=(const state &s)
{
	state snew = *this;
	return snew *= s;
}

state state::operator-() const
{
	state s = *this;
	for(auto iter = s.begin(); iter != s.end(); iter++)
		iter->second = -iter->second;
	return s;
}

state state::operator*(complex_type x) const
{
	state s = *this;
	return s *= x;
}

state linopt::operator*(complex_type x, const state &s)
{
	return s*x;
}

state &state::operator*=(complex_type x)
{
	for(auto iter = begin(); iter != end(); iter++)
		iter->second *= x;
	return *this;
}

state state::operator/(complex_type x) const
{
	state s = *this;
	return s /= x;
}

state &state::operator/=(complex_type x)
{
	for(auto iter = begin(); iter != end(); iter++)
		iter->second /= x;
	return *this;
}

real_type state::norm() const
{
	real_type n = 0.;
	for(auto iter = begin(); iter != end(); iter++)
		n += std::norm(iter->second);
	return sqrt(n);
}

state &state::normalize()
{
	return *this /= norm();
}

complex_type state::dot(const state &s) const
{
	complex_type z = 0.;
	auto a = this->begin();
	auto b = s.begin();
	while(a != this->end() && b != s.end())
	{
		if(a->first < b->first)
		{
			a++;
		}
		else if(a->first == b->first)
		{
			z += conj(a->second) * (b->second);
			a++;
			b++;
		}
		else
		{
			b++;
		}
	}
	return z;
}

complex_type linopt::dot(const state &a, const state &b)
{
	return a.dot(b);
}

state state::postselect(const fock &ancilla) const
{
	state s;
	auto asize = ancilla.size();
	for(auto iter = begin(); iter != end(); iter++)
	{
		const fock &f = iter->first;
		if(std::equal(ancilla.rbegin(), ancilla.rend(), f.rbegin()))
			s.emplace_hint(s.end(), fock(f.begin(), f.end() - asize), iter->second);
	}
	return s;
}

basis state::get_basis() const
{
	basis b;
	for(auto iter = begin(); iter != end(); iter++)
		b.insert(b.end(), iter->first);
	return b;
}

template<typename T>
std::ostream& print_array(std::ostream &stream, const T &a,
						  const char *b1 = "{",
						  const char *delim = ", ",
						  const char *b2 = "}")
{
	stream << b1;
	for(auto iter = a.begin(); iter != --a.end(); iter++)
		stream << *iter << delim;
	stream << *(--a.end()) << b2;
	return stream;
}

std::ostream& operator<<(std::ostream &stream, const linopt::state_element &e)
{
	stream << e.first << " = " << e.second;
	return stream;
}

std::ostream& operator<<(std::ostream &stream, const fock &f)
{
	return print_array(stream, f);
}

std::ostream& operator<<(std::ostream &stream, const basis &b)
{
	return print_array(stream, b, "{", ",\n", "}");
}

std::ostream& operator<<(std::ostream &stream, const state &s)
{
	return print_array(stream, s, "{", ",\n", "}");
}

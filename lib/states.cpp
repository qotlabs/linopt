#include "misc.h"
#include "states.h"

#include <algorithm>

using namespace linopt;

int fock::total() const
{
	int tot = 0;
	for(auto &n: *this)
		tot += n;
	return tot;
}

static constexpr int factorial_precomputed[] = {
	1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600};

static inline real_type factorial(int n)
{
	return (n <= 12) ? factorial_precomputed[n] : std::tgamma(n+1);
}

real_type fock::prod_fact() const
{
	real_type p = 1.;
	for(auto &n: *this)
		p *= factorial(n);
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

basis::basis(int nphot, int modes):
	basis()
{
	generate_basis(nphot, modes);
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
	for(auto &f1: *this)
		for(auto &f2: b)
			newb.insert(newb.end(), f1 * f2);
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
	for(auto &fp: *this)
	{
		if(std::equal(ancilla.begin(), ancilla.end(), fp.begin()))
			b.insert(b.end(), fock(fp.begin() + ancilla.size(), fp.end()));
	}
	return b;
}

state basis::apply_func(const basis_func &f) const
{
	state s;
	for(auto &elem: *this)
		s.insert(s.end(), state::element(elem, f(elem)));
	return s;
}

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
			snew.insert(snew.end(), state::element(a->first, a->second + b->second));
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
	for(auto &elem: s)
		(*this)[elem.first] += elem.second;
	return *this;
}

state &state::operator-=(const state &s)
{
	for(auto &elem: s)
		(*this)[elem.first] -= elem.second;
	return *this;
}

state state::operator*(const state &s) const
{
	state snew;
	for(auto &a: *this)
		for(auto &b: s)
			snew.insert(snew.end(), state::element(a.first * b.first, a.second * b.second));
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
	for(auto &elem: s)
		elem.second = -elem.second;
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
	for(auto &elem: *this)
		elem.second *= x;
	return *this;
}

state state::operator/(complex_type x) const
{
	state s = *this;
	return s /= x;
}

state &state::operator/=(complex_type x)
{
	for(auto &elem: *this)
		elem.second /= x;
	return *this;
}

real_type state::norm() const
{
	real_type n = 0.;
	for(auto &elem: *this)
		n += std::norm(elem.second);
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
	bool found = false;
	for(auto &elem: *this)
	{
		const fock &f = elem.first;
		const complex_type &amp = elem.second;
		if(std::equal(ancilla.begin(), ancilla.end(), f.begin()))
		{
			s.emplace_hint(s.end(), fock(f.begin() + asize, f.end()), amp);
			found = true;
		}
		else if(found)
		{
			break;
		}
	}
	return s;
}

std::map<fock, state> state::postselect(int modes) const
{
	std::map<fock, state> res;
	const fock &f = (*this->begin()).first;
	fock anc(f.begin(), f.begin() + modes);
	state *s = &res[anc];
	for(auto &elem: *this)
	{
		const fock &f = elem.first;
		const complex_type &amp = elem.second;
		if(!std::equal(anc.begin(), anc.end(), f.begin()))
		{
			anc.assign(f.begin(), f.begin() + modes);
			s = &res[anc];
		}
		s->emplace_hint(s->end(), fock(f.begin() + modes, f.end()), amp);
	}
	return res;
}

basis state::get_basis() const
{
	basis b;
	for(auto &elem: *this)
		b.insert(b.end(), elem.first);
	return b;
}

std::ostream& operator<<(std::ostream &stream, const linopt::state::element &e)
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

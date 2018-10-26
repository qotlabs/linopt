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

/** @defgroup states States
 * @brief Linear optics states.
 */

#include "misc.h"
#include "states.h"
#include "exceptions.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <memory>

using namespace linopt;

/**
 * @brief Returns the total number of photons in all modes.
 */
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

/**
 * @brief Returns a product of factorials of occupation numbers.
 */
real_type fock::prod_fact() const
{
	real_type p = 1.;
	for(auto &n: *this)
		p *= factorial(n);
	return p;
}

/**
 * @brief Returns a tensor product of `*this` and `f`.
 */
fock fock::operator*(const fock &f) const
{
	fock newf = *this;
	return newf *= f;
}

/**
 * @brief Effectively equivalent to `*this = (*this) * f`.
 */
fock &fock::operator*=(const fock &f)
{
	insert(end(), f.begin(), f.end());
	return *this;
}

/**
 * @brief Returns a sum of two Fock states (elementwise addition of
 * corresponding occupation numbers).
 *
 * @throw
 * If `*this` and `f` have different sizes then `wrong_size` is thrown.
 */
fock fock::operator+(const fock &f) const
{
	fock newf = *this;
	return newf += f;
}

/**
 * @brief Effectively equivalent to `*this = (*this) + f`.
 */
fock &fock::operator+=(const fock &f)
{
	if(this->size() != f.size())
		throw wrong_size(ERROR_MSG("Sizes of two Fock states should be equal. "
				"Currently they are " + std::to_string(this->size()) + " and " +
				std::to_string(f.size()) + "."));
	auto iter = this->begin();
	for(auto &n: f)
	{
		*iter += n;
		iter++;
	}
	return *this;
}

/**
 * @brief Constructs a basis of all possible Fock states with `modes` modes and
 * containing `nphot` photons.
 */
basis::basis(int nphot, int modes):
	basis()
{
	generate_basis(nphot, modes);
}

/**
 * @brief Returns a basis which is a union of Fock states from both `*this` and
 * `b`.
 */
basis basis::operator+(const basis &b) const
{
	basis newb = *this;
	return newb += b;
}

/**
 * @brief Effectively equivalent to `*this = (*this) + b`.
 */
basis &basis::operator+=(const basis &b)
{
	insert(b.begin(), b.end());
	return *this;
}

/**
 * @brief Calculates a tensor product of two bases.
 *
 * Returns a basis consisting of all possible elementwise tensor
 * products of elements of `*this` and `b`.
 */
basis basis::operator*(const basis &b) const
{
	basis newb;
	for(auto &f1: *this)
		for(auto &f2: b)
			newb.insert(newb.end(), f1 * f2);
	return newb;
}

/**
 * @brief Effectively equivalent to `*this = (*this) * b`.
 */
basis &basis::operator*=(const basis &b)
{
	*this = (*this) * b;
	return *this;
}

/**
 * @brief Generates a basis of all possible Fock states with `modes` modes and
 * containing `nphot` photons.
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
 * `basis::clear` first.
 */
basis &basis::generate_basis(const int nphot, const int modes, const fock &head)
{
	if(head.size() == modes)
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

/**
 * @brief Returns a postselected basis after observing ancilla.
 *
 * @param[in] ancilla -- ancilla's Fock state for postselection.
 * @return A basis after postselection.
 *
 * Postselection of a basis @f$ B @f$ with the ancilla state
 * @f$ | \mathrm{anc} \rangle = | a_0, a_2, \dots, a_A \rangle @f$
 * picks only Fock states from @f$ B @f$ with the heading (first @f$ A @f$
 * occupation numbers) which equals to @f$ | \mathrm{anc} \rangle @f$.
 * The headings do not get to the new constructed basis.
 * Note that currently there is no possibility to specify a position of ancilla
 * modes.
 *
 * For example, postselection of a basis
 * @f[ \{ | 0000001 \rangle, | 1234567 \rangle, | 1239999 \rangle \} @f]
 * with the ancilla @f$ | 123 \rangle @f$ results in the basis
 * @f[ \{ | 4567 \rangle, | 9999 \rangle \}. @f]
 *
 */
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

/**
 * @brief Constructs a state from the basis using function `f`.
 *
 * @param[in] f -- function to apply.
 * @return Constructed state.
 *
 * Applies a function `f` to all Fock states of `*this` to compute a
 * corresponding amplitude of a resulting state.
 *
 * @note
 * This method is automatically parallelized using OpenMP.
 */
template<class exec_policy>
state basis::apply_function(const fock_amp_function &f) const
{
	state s(*this);
	s.set_amplitudes<exec_policy>(f);
	return s;
}

template
state basis::apply_function<execution::seq>(const fock_amp_function &f) const;

template
state basis::apply_function<execution::par>(const fock_amp_function &f) const;

/**
 * @brief Adds two states, i.e., calculates their superposition.
 */
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

/**
 * @brief Effectively equivalent to `*this = (*this) + s`.
 */
state &state::operator+=(const state &s)
{
	for(auto &elem: s)
		(*this)[elem.first] += elem.second;
	return *this;
}

/**
 * @brief Effectively equivalent to `*this = (*this) - s`.
 */
state &state::operator-=(const state &s)
{
	for(auto &elem: s)
		(*this)[elem.first] -= elem.second;
	return *this;
}

/**
 * @brief Returns a tensor product of two states.
 */
state state::operator*(const state &s) const
{
	state snew;
	for(auto &a: *this)
		for(auto &b: s)
			snew.insert(snew.end(), state::element(a.first * b.first, a.second * b.second));
	return snew;
}

/**
 * @brief Effectively equivalent to `*this = (*this) * s`.
 */
state &state::operator*=(const state &s)
{
	state snew = *this;
	return snew *= s;
}

/**
 * @brief Negates amplitudes of the state.
 */
state state::operator-() const
{
	state s = *this;
	for(auto &elem: s)
		elem.second = -elem.second;
	return s;
}

/**
 * @brief Multiplies a state by a complex number.
 */
state state::operator*(complex_type x) const
{
	state s = *this;
	return s *= x;
}

/** @ingroup states
 * @brief Multiplies a state by a complex number.
 */
state linopt::operator*(complex_type x, const state &s)
{
	return s*x;
}

/**
 * @brief Effectively equivalent to `*this = (*this) * x`.
 */
state &state::operator*=(complex_type x)
{
	for(auto &elem: *this)
		elem.second *= x;
	return *this;
}

/**
  * @brief Divides a state by a complex number.
  */
state state::operator/(complex_type x) const
{
	state s = *this;
	return s /= x;
}

/**
 * @brief Effectively equivalent to `*this = (*this) / x`.
 */
state &state::operator/=(complex_type x)
{
	for(auto &elem: *this)
		elem.second /= x;
	return *this;
}

/**
 * @brief Returns norm of the state.
 */
real_type state::norm() const
{
	real_type n = 0.;
	for(auto &elem: *this)
		n += std::norm(elem.second);
	return sqrt(n);
}

/**
 * @brief Normalizes the state to have unit norm.
 */
state &state::normalize()
{
	return *this /= norm();
}

/**
 * @brief Calculates a dot (scalar) product.
 */
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

/** @ingroup states
 * @brief Calculates a dot (scalar) product.
 */
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

std::map<fock, state> state::postselect(const basis &b) const
{
	if(b.size() == 0)
		return {{fock(), *this}};
	std::map<fock, state> res;
	auto si = this->begin();
	auto bi = b.begin();
	state *res_bi = &res[*bi];
	while(si != this->end() && bi != b.end())
	{
		const fock &f = si->first;
		const complex_type &amp = si->second;
		const fock anc(f.begin(), f.begin() + bi->size());
		if(*bi < anc)
		{
			bi++;
			if(bi == b.end())
				break;
			res_bi = &res[*bi];
		}
		else if(*bi == anc)
		{
			res_bi->emplace_hint(res_bi->end(), fock(f.begin() + bi->size(), f.end()), amp);
			si++;
		}
		else
		{
			si++;
		}
	}
	bi++;
	for(; bi != b.end(); bi++)
		res[*bi] = state();
	return res;
}

basis state::get_basis() const
{
	basis b;
	for(auto &elem: *this)
		b.insert(b.end(), elem.first);
	return b;
}

void state::set_basis(const basis &b)
{
	clear();
	for(auto &f: b)
		emplace_hint(end(), f, 0.);
}

std::vector<state::value_type> state::get_amplitudes() const
{
	std::vector<state::value_type> amps;
	amps.reserve(size());
	for(auto &elem: *this)
		amps.push_back(elem.second);
	return amps;
}

template<>
void state::set_amplitudes<execution::seq>(const std::vector<complex_type> &amps)
{
	using std::to_string;
	if(static_cast<int>(amps.size()) != size())
		throw wrong_size(ERROR_MSG("The size of 'amps' (which is " +
			to_string(amps.size()) + ") should be equal to the state size (which is " +
			to_string(size()) + ")."));
	auto amps_iter = amps.begin();
	for(auto &elem: *this)
	{
		elem.second = *amps_iter;
		amps_iter++;
	}
}

template<>
void state::set_amplitudes<execution::par>(const std::vector<complex_type> &amps)
{
	throw not_implemented(ERROR_MSG("Parallel execution is not supported."));
}

template<>
void state::set_amplitudes<execution::seq>(const fock_amp_function &f)
{
	for(auto &elem: *this)
		elem.second = f(elem.first);
}

template<>
void state::set_amplitudes<execution::par>(const fock_amp_function &f)
{
	const int tnum = omp_get_max_threads();
	#pragma omp parallel num_threads(tnum)
	{
		const int tid = omp_get_thread_num();
		int cnt = 0;
		for(auto iter = begin(); iter != end(); iter++, cnt++)
			if(cnt % tnum == tid)
				iter->second = f(iter->first);
	}
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

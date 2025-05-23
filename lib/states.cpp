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

/** @defgroup states States
 * @brief Linear optics states.
 */

#include "states.h"
#include "misc.h"
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
int Fock::total() const
{
	int tot = 0;
	for(auto &n: *this)
		tot += n;
	return tot;
}

static constexpr int factorialPrecomputed[] = {
	1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600};

static inline Real factorial(int n)
{
	return (n <= 12) ? factorialPrecomputed[n] : std::tgamma(n+1);
}

/**
 * @brief Returns a product of factorials of occupation numbers.
 */
Real Fock::prodFact() const
{
	Real p = 1.;
	for(auto &n: *this)
		p *= factorial(n);
	return p;
}

/**
 * @brief Returns a tensor product of `*this` and `f`.
 */
Fock Fock::operator*(const Fock &f) const
{
	Fock newf = *this;
	return newf *= f;
}

/**
 * @brief Effectively equivalent to `*this = (*this) * f`.
 */
Fock &Fock::operator*=(const Fock &f)
{
	insert(end(), f.begin(), f.end());
	return *this;
}

/**
 * @brief Returns a sum of two Fock states (elementwise addition of
 * corresponding occupation numbers).
 *
 * @throw
 * If `*this` and `f` have different sizes then `WrongSize` is thrown.
 */
Fock Fock::operator+(const Fock &f) const
{
	Fock newf = *this;
	return newf += f;
}

/**
 * @brief Effectively equivalent to `*this = (*this) + f`.
 */
Fock &Fock::operator+=(const Fock &f)
{
	if(this->size() != f.size())
		throw WrongSize(ERROR_MSG("Sizes of two Fock states should be equal. "
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
Basis::Basis(int nphot, int modes):
	Basis()
{
	generateBasis(nphot, modes);
}

/**
 * @brief Returns a basis which is a union of Fock states from both `*this` and
 * `b`.
 */
Basis Basis::operator+(const Basis &b) const
{
	Basis newb = *this;
	return newb += b;
}

/**
 * @brief Effectively equivalent to `*this = (*this) + b`.
 */
Basis &Basis::operator+=(const Basis &b)
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
Basis Basis::operator*(const Basis &b) const
{
	Basis newb;
	for(auto &f1: *this)
		for(auto &f2: b)
			newb.insert(newb.end(), f1 * f2);
	return newb;
}

/**
 * @brief Effectively equivalent to `*this = (*this) * b`.
 */
Basis &Basis::operator*=(const Basis &b)
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
Basis &Basis::generateBasis(const int nphot, const int modes, const Fock &head)
{
	if(head.size() == modes)
	{
		if(nphot == 0)
			this->insert(this->end(), head);
		return *this;
	}
	for(int i = 0; i <= nphot; i++)
	{
		Fock f(head);
		f.pushBack(i);
		generateBasis(nphot - i, modes, f);
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
Basis Basis::postselect(const Fock &ancilla) const
{
	Basis b;
	for(auto &fp: *this)
	{
		if(std::equal(ancilla.begin(), ancilla.end(), fp.begin()))
			b.insert(b.end(), Fock(fp.begin() + ancilla.size(), fp.end()));
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
 */
template<typename ExecPolicy>
State Basis::applyFunction(const FockAmpFunction &f) const
{
	State s(*this);
	s.setAmplitudes<ExecPolicy>(f);
	return s;
}

template State Basis::applyFunction<execution::Seq>(const FockAmpFunction &f) const;
template State Basis::applyFunction<execution::Par>(const FockAmpFunction &f) const;

/**
 * @brief Adds two states, i.e., calculates their superposition.
 */
State State::operator+(const State &s) const
{
	State snew;
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
			snew.insert(snew.end(), State::Element(a->first, a->second + b->second));
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
State &State::operator+=(const State &s)
{
	for(auto &elem: s)
		(*this)[elem.first] += elem.second;
	return *this;
}

/**
 * @brief Effectively equivalent to `*this = (*this) - s`.
 */
State &State::operator-=(const State &s)
{
	for(auto &elem: s)
		(*this)[elem.first] -= elem.second;
	return *this;
}

/**
 * @brief Returns a tensor product of two states.
 */
State State::operator*(const State &s) const
{
	State snew;
	for(auto &a: *this)
		for(auto &b: s)
			snew.insert(snew.end(), State::Element(a.first * b.first, a.second * b.second));
	return snew;
}

/**
 * @brief Effectively equivalent to `*this = (*this) * s`.
 */
State &State::operator*=(const State &s)
{
	return *this = (*this) * s;
}

/**
 * @brief Negates amplitudes of the state.
 */
State State::operator-() const
{
	State s = *this;
	for(auto &elem: s)
		elem.second = -elem.second;
	return s;
}

/**
 * @brief Multiplies a state by a complex number.
 */
State State::operator*(Complex x) const
{
	State s = *this;
	return s *= x;
}

/** @ingroup states
 * @brief Multiplies a state by a complex number.
 */
State linopt::operator*(Complex x, const State &s)
{
	return s*x;
}

/**
 * @brief Effectively equivalent to `*this = (*this) * x`.
 */
State &State::operator*=(Complex x)
{
	for(auto &elem: *this)
		elem.second *= x;
	return *this;
}

/**
  * @brief Divides a state by a complex number.
  */
State State::operator/(Complex x) const
{
	State s = *this;
	return s /= x;
}

/**
 * @brief Effectively equivalent to `*this = (*this) / x`.
 */
State &State::operator/=(Complex x)
{
	for(auto &elem: *this)
		elem.second /= x;
	return *this;
}

/**
 * @brief Returns norm of the state.
 */
Real State::norm() const
{
	Real n = 0.;
	for(auto &elem: *this)
		n += std::norm(elem.second);
	return sqrt(n);
}

/**
 * @brief Normalizes the state to have unit norm.
 */
State &State::normalize()
{
	return *this /= norm();
}

/**
 * @brief Calculates a dot (scalar) product.
 */
Complex State::dot(const State &s) const
{
	Complex z = 0.;
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
Complex linopt::dot(const State &a, const State &b)
{
	return a.dot(b);
}

State State::postselect(const Fock &ancilla) const
{
	State s;
	auto asize = ancilla.size();
	bool found = false;
	for(auto &elem: *this)
	{
		const Fock &f = elem.first;
		const Complex &amp = elem.second;
		if(std::equal(ancilla.begin(), ancilla.end(), f.begin()))
		{
			s.emplace_hint(s.end(), Fock(f.begin() + asize, f.end()), amp);
			found = true;
		}
		else if(found)
		{
			break;
		}
	}
	return s;
}

std::map<Fock, State> State::postselect(int modes) const
{
	std::map<Fock, State> res;
	const Fock &f = (*this->begin()).first;
	Fock anc(f.begin(), f.begin() + modes);
	State *s = &res[anc];
	for(auto &elem: *this)
	{
		const Fock &f = elem.first;
		const Complex &amp = elem.second;
		if(!std::equal(anc.begin(), anc.end(), f.begin()))
		{
			anc.assign(f.begin(), f.begin() + modes);
			s = &res[anc];
		}
		s->emplace_hint(s->end(), Fock(f.begin() + modes, f.end()), amp);
	}
	return res;
}

std::map<Fock, State> State::postselect(const Basis &b) const
{
	if(b.size() == 0)
		return {{Fock(), *this}};
	std::map<Fock, State> res;
	auto si = this->begin();
	auto bi = b.begin();
	State *resBi = &res[*bi];
	while(si != this->end() && bi != b.end())
	{
		const Fock &f = si->first;
		const Complex &amp = si->second;
		const Fock anc(f.begin(), f.begin() + bi->size());
		if(*bi < anc)
		{
			bi++;
			if(bi == b.end())
				break;
			resBi = &res[*bi];
		}
		else if(*bi == anc)
		{
			resBi->emplace_hint(resBi->end(), Fock(f.begin() + bi->size(), f.end()), amp);
			si++;
		}
		else
		{
			si++;
		}
	}
	bi++;
	for(; bi != b.end(); bi++)
		res[*bi] = State();
	return res;
}

Basis State::getBasis() const
{
	Basis b;
	for(auto &elem: *this)
		b.insert(b.end(), elem.first);
	return b;
}

void State::setBasis(const Basis &b)
{
	clear();
	for(auto &f: b)
		emplace_hint(end(), f, 0.);
}

std::vector<State::Value> State::getAmplitudes() const
{
	std::vector<State::Value> amps;
	amps.reserve(size());
	for(auto &elem: *this)
		amps.push_back(elem.second);
	return amps;
}

template<>
void State::setAmplitudes<execution::Seq>(const std::vector<Complex> &amps)
{
	using std::to_string;
	if(static_cast<int>(amps.size()) != size())
		throw WrongSize(ERROR_MSG("The size of 'amps' (which is " +
			to_string(amps.size()) + ") should be equal to the state size (which is " +
			to_string(size()) + ")."));
	auto ampsIter = amps.begin();
	for(auto &elem: *this)
	{
		elem.second = *ampsIter;
		ampsIter++;
	}
}

template<>
void State::setAmplitudes<execution::Par>(const std::vector<Complex> &)
{
	throw NotSupported(ERROR_MSG("Parallel execution is not supported."));
}

template<>
void State::setAmplitudes<execution::Seq>(const FockAmpFunction &f)
{
	for(auto &elem: *this)
		elem.second = f(elem.first);
}

#ifdef _OPENMP
template<>
void State::setAmplitudes<execution::Par>(const FockAmpFunction &f)
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

#else
template<>
void state::setAmplitudes<execution::par>(const FockAmpFunction &)
{
	throw NotSupported(ERROR_MSG("Parallel execution is not supported when "
								 "the library is compiled without OpenMP."));
}
#endif // _OPENMP

std::ostream& operator<<(std::ostream &stream, const linopt::State::Element &e)
{
	stream << e.first << " = " << e.second;
	return stream;
}

std::ostream& operator<<(std::ostream &stream, const Fock &f)
{
	return printArray(stream, f);
}

std::ostream& operator<<(std::ostream &stream, const Basis &b)
{
	return printArray(stream, b, "{", ",\n", "}");
}

std::ostream& operator<<(std::ostream &stream, const State &s)
{
	return printArray(stream, s, "{", ",\n", "}");
}

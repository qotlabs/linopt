#include "linopt.h"

#include <algorithm>

using namespace linopt;

basis::basis():
    std::set<fock>() {}

basis::basis(const basis &b):
    std::set<fock>(b) {}

basis::basis(std::initializer_list<fock> il):
    std::set<fock>(il) {}

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

basis &basis::generate_basis(int nphot, int modes)
{
    return *this;
}

basis basis::postselect(const fock &ancilla) const
{
    basis b;
    for(auto fp = begin(); fp != end(); fp++)
    {
        if(std::equal(ancilla.rbegin(), ancilla.rend(), fp->rbegin()))
            b.insert(b.end(),
                     fock(fp->begin(), fp->begin() + fp->size() - ancilla.size()));
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

state state::operator*(complex_type x) const
{
    state s = *this;
    return s *= x;
}

state &state::operator*=(complex_type x)
{
    for(auto iter = begin(); iter != end(); iter++)
        iter->second *= x;
    return *this;
}

state operator/(complex_type x) const
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
    return 0.;
}

state state::postselect(const fock &ancilla) const
{
    state s;
    for(auto iter = begin(); iter != end(); iter++)
    {
        state_element e = *iter;
        fock f = e.first;
        complex_type val = e.second;
        if(std::equal(ancilla.rbegin(), ancilla.rend(), f.rbegin()))
        {
            f = fock(f.begin(), f.begin() + f.size() - ancilla.size());
            s.insert(s.end(), state_element(f, val));
        }
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
std::ostream& print_array(std::ostream &stream, const T &a)
{
    stream << "{";
    for(auto iter = a.begin(); iter != --a.end(); iter++)
        stream << *iter << ", ";
    stream << *(--a.end()) << "}";
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
    return print_array(stream, b);
}

std::ostream& operator<<(std::ostream &stream, const state &s)
{
    return print_array(stream, s);
}

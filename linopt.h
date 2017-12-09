#ifndef LINOPT_H
#define LINOPT_H

#include <vector>
#include <set>
#include <map>
#include <ostream>
#include <initializer_list>
#include <functional>
#include <complex>

namespace linopt
{

class basis;
class state;

typedef double real_type;
typedef std::complex<real_type> complex_type;
typedef std::vector<int> fock;
typedef std::pair<fock, complex_type> state_element;
typedef std::function<complex_type(fock)> basis_func;

class basis : public std::set<fock>
{
public:
    basis();
    basis(const basis &b);
    basis(std::initializer_list<fock> il);
    explicit basis(int nphot, int modes);
    basis operator+(const basis &b) const;
    basis &operator+=(const basis &b);
    basis &generate_basis(int nphot, int modes);    // TODO
    basis postselect(const fock &ancilla) const;
    state apply_func(const basis_func &f) const;
};

class state : public std::map<fock, complex_type>
{
public:
    state();
    state(const state &s);
    state operator+(const state &s) const;  // TODO
    state &operator+=(const state &s);      // TODO
    state operator-(const state &s) const;  // TODO
    state &operator-=(const state &s);      // TODO
    state operator*(complex_type x) const;
    state &operator*=(complex_type x);
    state operator/(complex_type x) const;
    state &operator/=(complex_type x);
    real_type norm() const;
    state &normalize();
    complex_type dot(const state &s) const; // TODO
    state postselect(const fock &ancilla) const;
    basis get_basis() const;
};

}

std::ostream &operator<<(std::ostream &stream, const linopt::fock &f);
std::ostream &operator<<(std::ostream &stream, const linopt::basis &b);
std::ostream &operator<<(std::ostream &stream, const linopt::state_element &e);
std::ostream &operator<<(std::ostream &stream, const linopt::state &s);

#endif // LINOPT_H

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

#include <iterator>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include "../lib/linopt.h"
#include "../lib/misc.h"

using namespace linopt;
namespace py = pybind11;

// __repr__ and __str__ methods for class "fock"
std::string fock_str(const fock &f)
{
	std::stringstream ss;
	print_array(ss, f, "[", ",", "]");
	return ss.str();
}

std::string fock_repr(const fock &f)
{
	return "fock(" + fock_str(f) + ")";
}

// __repr__ and __str__ methods for class "basis"
std::string basis_str(const basis &b)
{
	std::stringstream ss;
	if(!b.empty())
	{
		ss << "{";
		for(auto iter = b.begin(); iter != --b.end(); ++iter)
			print_array(ss, *iter, "(", ",", ")") << ", ";
		print_array(ss, *(--b.end()), "(", ",", ")");
		ss << "}";
	}
	else
	{
		ss << "{}";
	}
	return ss.str();
}

std::string basis_repr(const basis &b)
{
	return "basis(" + basis_str(b) + ")";
}

// __repr__ and __str__ methods for class "state"
std::string state_str(const state &s)
{
	std::stringstream ss;
	if(!s.empty())
	{
		ss << "{";
		for(auto iter = s.begin(); iter != --s.end(); ++iter)
		{
			print_array(ss, iter->first, "(", ",", ")") << ": ";
			print_complex(ss, iter->second);
			ss << ", ";
		}
		auto iter = --s.end();
		print_array(ss, iter->first, "(", ",", ")") << ": ";
		print_complex(ss, iter->second);
		ss << "}";
	}
	else
	{
		ss << "{}";
	}
	return ss.str();
}

std::string state_repr(const state &s)
{
	return "state(" + state_str(s) + ")";
}

template<typename Container, typename Key>
void erase_key(Container &c, const Key &k)
{
	auto iter = c.find(k);
	if(iter == c.end())
		throw py::key_error();
	else
		c.erase(k);
}

static const auto docstr =
	"Linear optics circuit calculator.";

#if PYBIND11_VERSION_MAJOR <= 2 && PYBIND11_VERSION_MINOR < 2
PYBIND11_PLUGIN(pylinopt)
{
	py::module m("pylinopt", docstr);
#else
PYBIND11_MODULE(pylinopt, m)
{
	m.doc() = docstr;
#endif

	// Matrix functions
	m
	.def("permanent", &permanent,
		 "Calculates permanent of a square nxn matrix M.\n"
		 "Internally the function uses Glynn formula with Gray code summation technique.\n"
		 "Complexity is O(n*2^n).",
		 py::arg("M"))

	.def("is_column_unitary", &is_column_unitary,
		 "Checks whether a matrix M is column unitary within accuracy eps.",
		 py::arg("M"), py::arg("eps") = default_epsilon)

	.def("is_row_unitary", &is_row_unitary,
		 "Checks whether a matrix M is row unitary within accuracy eps.",
		 py::arg("M"), py::arg("eps") = default_epsilon)

	.def("is_unitary", &is_unitary,
		 "Checks whether a square matrix M is unitary within accuracy eps.",
		 py::arg("M"), py::arg("eps") = default_epsilon)

	.def("hurwitz", (matrix_type (*)(const point &x)) &hurwitz_parametrization,
		 "Hurwitz parametrization of a unitary matrix.",
		 py::arg("x"))

	.def("exp_hermite", (matrix_type (*)(const point &x)) &exp_hermite_parametrization,
		 "Exponential Hermitian parametrization of a unitary matrix.\n"
		 "The resulting matrix is, U = exp(j*H(x)).",
		 py::arg("x"))
	;

	// Fock
	py::class_<fock>(m, "fock")
		.def(py::init<>())
		.def(py::init<const fock &>())
		.def(py::init<const fock::vector_class &>())

		.def("__str__", &fock_str)
		.def("__repr__", &fock_repr)

		.def(py::self < py::self,
			 "Compares two Fock states in lexicographic order. Returns True if "
			 "'self' is less than 'f'.",
			 py::arg("f"))
		.def(py::self <= py::self,
			 "Compares two Fock states in lexicographic order. Returns True if "
			 "'self' is less or equal than 'f'.",
			 py::arg("f"))
		.def(py::self == py::self,
			 "Tests whether two Fock states are equal.",
			 py::arg("f"))
		.def(py::self != py::self,
			 "Tests whether two Fock state differ.",
			 py::arg("f"))
		.def(py::self >= py::self,
			 "Compares two Fock states in lexicographic order. Returns True if "
			 "'self' is greater or equal than 'f'.",
			 py::arg("f"))
		.def(py::self > py::self,
			 "Compares two Fock states in lexicographic order. Returns True if "
			 "'self' is greater than 'f'.",
			 py::arg("f"))

		.def("__len__", [](const fock &f) { return f.size(); },
			 "Returns number of modes of the Fock state.")

		.def("__getitem__", [](const fock &f, int i) {
				if( !(0 <= i && i < f.size()) )
					throw py::index_error("Index " + std::to_string(i) + " is out of range.");
				return f[i];
			 },
			 "Returns the occupation number of the i-th mode.",
			 py::arg("i"))

		.def("__getitem__", [](const fock &f, py::slice slice) {
				size_t start, stop, step, len;
				if(!slice.compute(f.size(), &start, &stop, &step, &len))
					throw py::error_already_set();
				fock *res = new fock(len);
				for(size_t i = 0; i < len; i++)
				{
					(*res)[i] = f[start];
					start += step;
				}
				return res;
			 })

		.def("__setitem__", [](fock &f, int i, fock::value_type val) {
				if( !(0 <= i && i < f.size()) )
					throw py::index_error("Index " + std::to_string(i) + " is out of range.");
				return f[i] = val;
			 },
			 "Sets the occupation number of the i-th mode to 'val'.",
			 py::arg("i"), py::arg("val"))

		.def("__setitem__", [](fock &f, py::slice slice, const fock &val) {
				 size_t start, stop, step, len;
				 if(!slice.compute(f.size(), &start, &stop, &step, &len))
					 throw py::error_already_set();
				 if(len != static_cast<size_t>(val.size()))
					 throw std::runtime_error("Left and right hand size of slice "
											  "assignment have different sizes.");
				 for(size_t i = 0; i < len; i++)
				 {
					 f[start] = val[i];
					 start += step;
				 }
			 })

		.def("__delitem__", [](fock &f, int i) {
				 if( !(0 <= i && i < f.size()) )
					 throw py::index_error("Index " + std::to_string(i) + " is out of range.");
				 return f.erase(f.begin() + i);
			 },
			 "Removes i-th mode.",
			 py::arg("i"))

		.def("__iter__", [](fock &f) {
				 return py::make_iterator(f.begin(), f.end());
			 },
			 "Returns corresponding iterator object.",
			 py::keep_alive<0, 1>())

		.def("__reversed__", [](fock &f) {
				 return py::make_iterator(f.rbegin(), f.rend());
			 },
			 "Returns corresponding reverse iterator object.",
			 py::keep_alive<0, 1>())

		.def("append", (void (fock::*)(const fock::value_type &n)) &fock::push_back,
			  "Adds the mode with the occupation number 'n' to the end.",
			  py::arg("n"))

		.def("insert", [](fock &f, int i, fock::value_type n) {
				if( !(0 <= i && i <= f.size()) )
					throw py::index_error("Index " + std::to_string(i) + " is out of range.");
				else
					f.insert(f.begin() + i, n);
			 },
			 "Inserts the mode with occupation number 'n' at position 'i'.",
			 py::arg("i"), py::arg("n"))

		.def("pop", [](fock &f) {
				 auto val = f.back();
				 f.pop_back();
				 return val;
			 },
			 "Removes the last mode and returns its occupation number.")

		.def("pop", [](fock &f, int i) {
				 if( !(0 <= i && i < f.size()) )
					 throw py::index_error("Index " + std::to_string(i) + " is out of range.");
				 auto val = f[i];
				 f.erase(f.begin() + i);
				 return val;
			 },
			 "Removes i-th mode and returns its occupation number.",
			 py::arg("i"))

		.def("clear", &fock::clear,
			 "Completely clears the Fock state")

		.def("resize", [](fock &f, int n) { f.resize(n); },
			 "Sets the numbers of modes for the Fock state to 'n'.",
			 py::arg("n"))

		.def("total", &fock::total,
			 "Calculates the total number of photons in all modes.")

		.def("prod_fact", &fock::prod_fact,
			 "Calculates a product of factorials of occupation numbers.")

		.def(py::self * py::self,
			 "Calculates a tensor product of two Fock states.")
		.def(py::self *= py::self,
			 "Effectively equivalent to self = self * f.",
			 py::arg("f"))

		.def(py::self + py::self,
			 "Calculates a sum of two Fock states (elementwise addition of "
			 "corresponding occupation numbers).")
		.def(py::self += py::self,
			 "Effectively equivalent to self = self + f.",
			 py::arg("f"))

		.def("as_list", [](const fock &f) {
				 py::list l;
				 for (const auto &elem: f)
					 l.append(elem);
				 return l;
			 },
			 "Returns Python list object that represents the Fock state.")
	;

	// Basis
	py::class_<basis>(m, "basis")
		.def(py::init<>())
		.def(py::init<const basis &>())
		.def(py::init<const basis::set_class &>())
		.def(py::init<int, int>(),
			 py::arg("nphot"), py::arg("modes"))

		.def("__str__", &basis_str)
		.def("__repr__", &basis_repr)

		.def("__len__", [](const basis &b) { return b.size(); },
			 "Returns the number of Fock states in the basis.")

		.def("__delitem__", (void (*)(basis &, const fock &)) &erase_key,
			 "Removes the Fock 'f' from the basis.",
			 py::arg("f"))

		.def("__iter__", [](basis &b) {
				 return py::make_iterator(b.begin(), b.end());
			 },
			 "Returns corresponding iterator object.",
			 py::keep_alive<0, 1>())

		.def("__reversed__", [](basis &b) {
				 return py::make_iterator(b.rbegin(), b.rend());
			 },
			 "Returns corresponding reverse iterator object.",
			 py::keep_alive<0, 1>())

		.def("__contains__", [](const basis &b, const fock &f) {
				 return b.find(f) != b.end();
			 },
			 "Tests whether the Fock 'f' is in the basis.",
			 py::arg("f"))

		.def("add", [](basis &b, const fock &f) {
				 b.insert(f);
			 },
			 "Adds the Fock state 'f' to the basis.",
			 py::arg("f"))

		.def("remove", (void (*)(basis &, const fock &)) &erase_key,
			 "Removes the Fock state 'f' from the basis. Throws KeyError if 'f' "
			 "does not exist",
			 py::arg("f"))

		.def("discard", [](basis &b, const fock &f) {
				 b.erase(f);
			 },
			 "Removes the Fock state 'f' from the basis if it is present.",
			 py::arg("f"))

		.def("clear", &basis::clear,
			 "Remove all elements from the basis.")

		.def(py::self + py::self,
			 "Returns a basis which is a union of Fock states from both bases.")
		.def(py::self += py::self,
			 "Effectively equivalent to self = self + b.",
			 py::arg("b"))

		.def(py::self * py::self,
			 "Calculates a tensor product of two bases.\n"
			 "Returns a basis consisting of all possible elementwise tensor "
			 "products of elements of the bases.")
		.def(py::self *= py::self,
			 "Effectively equivalent to self = self * b.",
			 py::arg("b"))

		.def("generate_basis", &basis::generate_basis,
			 "Generates a basis of all possible Fock states with 'modes' modes "
			 "and containing 'nphot' photons.",
			 py::arg("nphot"), py::arg("modes"), py::arg("head") = fock())

		.def("postselect", &basis::postselect,
			 "Returns a postselected basis after observing ancilla 'anc'.\n"
			 "Ancilla is assumed to occupy the first modes.",
			 py::arg("anc"))

		.def("apply_function", &basis::apply_function,
			 "Applies a function 'func' to all Fock states of the basis to "
			 "compute a corresponding amplitude and returns the corresponding "
			 "state. The 'func' should take a fock object as an argument and "
			 "return a complex number representing its amplitude.",
			 py::arg("func"))

		.def("as_set", [](const basis &b) {
				 py::set s;
				 for(const auto &f: b)
					 s.add(f);
				 return s;
			 },
			 "Returns Python set object that represents the basis.")
	;

	// State
	py::class_<state>(m, "state")
		.def(py::init<>())
		.def(py::init<const state &>())
		.def(py::init<const fock &>())
		.def(py::init<const state::map_class &>())

		.def("__str__", &state_str)
		.def("__repr__", &state_repr)

		.def("__len__", [](const state &s) { return s.size(); },
			 "Returns number of Focks with a specified amplitude in the state.")

		.def("__getitem__", [](const state &s, const fock &f) {
				 auto iter = s.find(f);
				 if(iter == s.end())
					 throw py::key_error();
				 return iter->second;
			 },
			 "Returns an amplitude corresponding to the Fock 'f'.",
			 py::arg("f"))

		.def("__missing__", [](state &s) { return state::value_type(0.); })

		.def("__setitem__", [](state &s, const fock &f, state::value_type amp) {
				 if(amp == 0.)
					 s.erase(f);
				 else
					 s[f] = amp;
				 return amp;
			 },
			 "Sets amplitude of the Fock 'f' to 'amp'. If 'amp' = 0, then "
			 "deletes the Fock 'f' from the state.",
			 py::arg("f"), py::arg("amp"))

		.def("__delitem__", (void (*)(state &, const fock &)) &erase_key,
			 "Removes the Fock 'f' from the state. Equivalent to self[f] = 0.",
			 py::arg("f"))

		.def("__iter__", [](state &s) {
				 return py::make_key_iterator(s.begin(), s.end());
			 },
			 "Returns corresponding iterator object.",
			 py::keep_alive<0, 1>())

		.def("__reversed__", [](state &s) {
				 return py::make_key_iterator(s.rbegin(), s.rend());
			 },
			 "Returns corresponding reverse iterator object.",
			 py::keep_alive<0, 1>())

		.def("__contains__", [](const state &s, const fock &f) {
				 return s.find(f) != s.end();
			 },
			 "Tests whether the state has a specified amplitude corresponding "
			 "to the Fock 'f'.",
			 py::arg("f"))

		.def("clear", &state::clear,
			 "Completely clears the state.")

		.def(py::self + py::self,
			 "Adds two states, i.e., calculates their superposition.")
		.def(py::self += py::self,
			 "Effectively equivalent to self = self + s.",
			 py::arg("s"))

		.def(py::self - py::self,
			 "Subtracts two states.")
		.def(py::self -= py::self,
			 "Effectively equivalent to self = self - s.",
			 py::arg("s"))

		.def(-py::self,
			 "Negates amplitudes of the state.")

		.def(py::self * py::self,
			 "Returns a tensor product of two states.")
		.def(py::self *= py::self,
			 "Effectively equivalent to self = self * s.",
			 py::arg("s"))

		.def(py::self * complex_type(),
			 "Multiplies a state by a complex number.")
		.def(py::self *= complex_type(),
			 "Effectively equivalent to self = self * z.",
			 py::arg("z"))

		.def(py::self / complex_type(),
			 "Divides a state by a complex number.")
		.def(py::self /= complex_type(),
			 "Effectively equivalent to self = self / z.",
			 py::arg("z"))

		.def("norm", &state::norm,
			 "Returns norm of the state.")

		.def("normalize", &state::normalize,
			 "Normalizes the state to have unit norm.")

		.def("dot", &state::dot,
			 "Calculates a dot (scalar) product.")

		.def("postselect", (state (state::*)(const fock&) const) &state::postselect,
			 "Returns postselected state for ancilla 'anc'. The ancilla is "
			 "assumed to occupy the first modes.",
			 py::arg("anc"))

		.def("postselect", (std::map<fock, state> (state::*)(int) const) &state::postselect,
			 "Calculates postselected states for all possible ancillas that "
			 "occupy the first 'nmodes' modes. Returns a dictionary in the "
			 "following format:\n"
			 "{anc1: postselected_state1, anc2: postselected_state2, ...},\n"
			 "where anc1, anc2, ... have 'nmodes' modes.\n"
			 "This function is useful when postselection for all possible "
			 "ancillas is required. It is faster than calling postselect(anc) "
			 "in cycle for different 'anc'.",
			 py::arg("nmodes"))

		.def("postselect", (std::map<fock, state> (state::*)(const basis &) const) &state::postselect,
			 "Calculates postselected states for all ancillas specified by "
			 "'basis'. Returns a dictionary in the following format:\n"
			 "{anc1: postselected_state1, anc2: postselected_state2, ...}.\n"
			 "This function is useful when postselection for an array of "
			 "ancillas is required. It is faster than calling postselect(anc) "
			 "in cycle for different 'anc'.",
			 py::arg("basis"))

		.def("get_basis", &state::get_basis,
			 "Returns basis of the state.")

		.def("as_dict", [](const state &s) {
				 py::dict d;
				 for(const auto &elem: s)
				 d[py::make_tuple(elem.first)] = elem.second;
			 return d;
			 },
			 "Returns Python dict object that represents the state.")
	;

	m
	.def("__mul__", (state (*)(complex_type, const state &)) &operator*,
		 "Multiplies a state by a complex number.",
		 py::arg("z"), py::arg("s"))
	.def("dot", &dot,
		 "Calculates a dot (scalar) product of two states.",
		 py::arg("s1"), py::arg("s2"))
	;

	// Implicit conversions
	py::implicitly_convertible<fock::vector_class, fock>();
	py::implicitly_convertible<basis::set_class, basis>();
	py::implicitly_convertible<state::map_class, state>();

	// Circuit
	py::class_<circuit>(m, "circuit")
		.def(py::init<>())

		.def_property("input_state",
			 &circuit::get_input_state,
			 &circuit::set_input_state,
			 "Input Fock state.")

		.def_property("output_basis",
			 &circuit::get_output_basis,
			 &circuit::set_output_basis,
			 "A basis of Fock states used to compute the output state. It may "
			 "either contain the full Fock state basis of the system or be "
			 "composed only of the combination of the possible states in the "
			 "logical and ancilla subsystems. Such a combination eliminates "
			 "the necessity to compute the probabilities for the output states "
			 "which won't be present in the system. For example, if the total "
			 "number of photons in the system is 4, the ancilla subsystem is "
			 "set to have 2 photons, and the logical subsystem is also set to "
			 "have 2 photons, the states, where 3 or 4 photons populate either "
			 "ancilla or logical subsystem are irrelevant to the problem being "
			 "solved.")

		.def_property("unitary",
			 &circuit::get_unitary,
			 &circuit::set_unitary,
			 "Unitary matrix representing a transformation of creation "
			 "operators of photons in modes.")

		.def("output_state", &circuit::output_state,
			 "Calculates output state.")
	;

#if PYBIND11_VERSION_MAJOR <= 2 && PYBIND11_VERSION_MINOR < 2
	return m.ptr();
#endif
}

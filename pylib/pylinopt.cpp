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

#include <iterator>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "matrix.h"
#include "states.h"
#include "circuit.h"
#include "circuit_design.h"
#include "misc.h"

using namespace linopt;
namespace py = pybind11;

enum class pyexecution
{
	seq,
	par
};

// __repr__ and __str__ methods for class "fock"
std::string fockStr(const Fock &f)
{
	std::stringstream ss;
	printArray(ss, f, "[", ",", "]");
	return ss.str();
}

std::string fockRepr(const Fock &f)
{
	return "fock(" + fockStr(f) + ")";
}

// __repr__ and __str__ methods for class "basis"
std::string basisStr(const Basis &b)
{
	std::stringstream ss;
	if(!b.empty())
	{
		ss << "{";
		for(auto iter = b.begin(); iter != --b.end(); ++iter)
			printArray(ss, *iter, "(", ",", ")") << ", ";
		printArray(ss, *(--b.end()), "(", ",", ")");
		ss << "}";
	}
	else
	{
		ss << "{}";
	}
	return ss.str();
}

std::string basisRepr(const Basis &b)
{
	return "basis(" + basisStr(b) + ")";
}

// __repr__ and __str__ methods for class "state"
std::string stateStr(const State &s)
{
	std::stringstream ss;
	if(!s.empty())
	{
		ss << "{";
		for(auto iter = s.begin(); iter != --s.end(); ++iter)
		{
			printArray(ss, iter->first, "(", ",", ")") << ": ";
			printComplex(ss, iter->second);
			ss << ", ";
		}
		auto iter = --s.end();
		printArray(ss, iter->first, "(", ",", ")") << ": ";
		printComplex(ss, iter->second);
		ss << "}";
	}
	else
	{
		ss << "{}";
	}
	return ss.str();
}

std::string stateRepr(const State &s)
{
	return "state(" + stateStr(s) + ")";
}

template<typename Container, typename Key>
void eraseKey(Container &c, const Key &k)
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
	py::module m("linopt", docstr);
#else
PYBIND11_MODULE(linopt, m)
{
	m.doc() = docstr;
#endif
	m.attr("default_epsilon") = py::float_(defaultEpsilon);

	// Execution policies
	py::enum_<pyexecution>(m, "execution")
		.value("seq", pyexecution::seq)
		.value("par", pyexecution::par)
		.export_values();

	// Matrix functions
	m
	.def("permanent", &permanent,
		 "Calculates permanent of a square nxn matrix M.\n"
		 "Internally the function uses Glynn formula with Gray code summation\n"
		 "technique.\n"
		 "Complexity is O(n*2^n).",
		 py::arg("M"))

	.def("is_column_unitary", &isColumnUnitary,
		 "Checks whether a matrix M is column unitary within accuracy eps.",
		 py::arg("M"), py::arg("eps") = defaultEpsilon)

	.def("is_row_unitary", &isRowUnitary,
		 "Checks whether a matrix M is row unitary within accuracy eps.",
		 py::arg("M"), py::arg("eps") = defaultEpsilon)

	.def("is_unitary", &isUnitary,
		 "Checks whether a square matrix M is unitary within accuracy eps.",
		 py::arg("M"), py::arg("eps") = defaultEpsilon)

	.def("hurwitz", (Matrix (*)(const Point &x)) &hurwitzParametrization,
		 "Hurwitz parametrization of a unitary matrix.",
		 py::arg("x"))

	.def("exp_hermite", (Matrix (*)(const Point &x)) &expHermiteParametrization,
		 "Exponential Hermitian parametrization of a unitary matrix.\n"
		 "The resulting matrix is, U = exp(j*H(x)).",
		 py::arg("x"))
	;

	// Fock
	py::class_<Fock>(m, "Fock")
		.def(py::init<>())
		.def(py::init<const Fock &>())
		.def(py::init<const Fock::Vector &>())

		.def("__str__", &fockStr)
		.def("__repr__", &fockRepr)

		.def(py::self < py::self,
			 "Compares two Fock states in lexicographic order. Returns True if\n"
			 "'self' is less than 'f'.",
			 py::arg("f"))
		.def(py::self <= py::self,
			 "Compares two Fock states in lexicographic order. Returns True if\n"
			 "'self' is less or equal than 'f'.",
			 py::arg("f"))
		.def(py::self == py::self,
			 "Tests whether two Fock states are equal.",
			 py::arg("f"))
		.def(py::self != py::self,
			 "Tests whether two Fock state differ.",
			 py::arg("f"))
		.def(py::self >= py::self,
			 "Compares two Fock states in lexicographic order. Returns True if\n"
			 "'self' is greater or equal than 'f'.",
			 py::arg("f"))
		.def(py::self > py::self,
			 "Compares two Fock states in lexicographic order. Returns True if\n"
			 "'self' is greater than 'f'.",
			 py::arg("f"))

		.def("__len__", [](const Fock &f) { return f.size(); },
			 "Returns number of modes of the Fock state.")

		.def("__getitem__", [](const Fock &f, int i) {
				if( !(0 <= i && i < f.size()) )
					throw py::index_error("Index " + std::to_string(i) + " is out of range.");
				return f[i];
			 },
			 "Returns the occupation number of the i-th mode.",
			 py::arg("i"))

		.def("__getitem__", [](const Fock &f, py::slice slice) {
				size_t start, stop, step, len;
				if(!slice.compute(f.size(), &start, &stop, &step, &len))
					throw py::error_already_set();
				Fock *res = new Fock(len);
				for(size_t i = 0; i < len; i++)
				{
					(*res)[i] = f[start];
					start += step;
				}
				return res;
			 })

		.def("__setitem__", [](Fock &f, int i, Fock::Value val) {
				if( !(0 <= i && i < f.size()) )
					throw py::index_error("Index " + std::to_string(i) + " is out of range.");
				return f[i] = val;
			 },
			 "Sets the occupation number of the i-th mode to 'val'.",
			 py::arg("i"), py::arg("val"))

		.def("__setitem__", [](Fock &f, py::slice slice, const Fock &val) {
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

		.def("__delitem__", [](Fock &f, int i) {
				 if( !(0 <= i && i < f.size()) )
					 throw py::index_error("Index " + std::to_string(i) + " is out of range.");
				 return f.erase(f.begin() + i);
			 },
			 "Removes i-th mode.",
			 py::arg("i"))

		.def("__iter__", [](Fock &f) {
				 return py::make_iterator(f.begin(), f.end());
			 },
			 "Returns corresponding iterator object.",
			 py::keep_alive<0, 1>())

		.def("__reversed__", [](Fock &f) {
				 return py::make_iterator(f.rbegin(), f.rend());
			 },
			 "Returns corresponding reverse iterator object.",
			 py::keep_alive<0, 1>())

		.def("append", &Fock::pushBack,
			  "Adds the mode with the occupation number 'n' to the end.",
			  py::arg("n"))

		.def("insert", [](Fock &f, int i, Fock::Value n) {
				if( !(0 <= i && i <= f.size()) )
					throw py::index_error("Index " + std::to_string(i) + " is out of range.");
				else
					f.insert(f.begin() + i, n);
			 },
			 "Inserts the mode with occupation number 'n' at position 'i'.",
			 py::arg("i"), py::arg("n"))

		.def("pop", [](Fock &f) {
				 auto val = f.back();
				 f.popBack();
				 return val;
			 },
			 "Removes the last mode and returns its occupation number.")

		.def("pop", [](Fock &f, int i) {
				 if( !(0 <= i && i < f.size()) )
					 throw py::index_error("Index " + std::to_string(i) + " is out of range.");
				 auto val = f[i];
				 f.erase(f.begin() + i);
				 return val;
			 },
			 "Removes i-th mode and returns its occupation number.",
			 py::arg("i"))

		.def("clear", [](Fock &f) { f.clear(); },
			 "Completely clears the Fock state")

		.def("resize", [](Fock &f, int n) { f.resize(n); },
			 "Sets the numbers of modes for the Fock state to 'n'.",
			 py::arg("n"))

		.def("total", &Fock::total,
			 "Calculates the total number of photons in all modes.")

		.def("prod_fact", &Fock::prodFact,
			 "Calculates a product of factorials of occupation numbers.")

		.def(py::self * py::self,
			 "Calculates a tensor product of two Fock states.")
		.def(py::self *= py::self,
			 "Effectively equivalent to self = self * f.",
			 py::arg("f"))

		.def(py::self + py::self,
			 "Calculates a sum of two Fock states (elementwise addition of\n"
			 "corresponding occupation numbers).")
		.def(py::self += py::self,
			 "Effectively equivalent to self = self + f.",
			 py::arg("f"))

		.def("as_list", [](const Fock &f) {
				 py::list l;
				 for (const auto &elem: f)
					 l.append(elem);
				 return l;
			 },
			 "Returns Python list object that represents the Fock state.")
	;

	// Basis
	py::class_<Basis>(m, "Basis")
		.def(py::init<>())
		.def(py::init<const Basis &>())
		.def(py::init<const Basis::Set &>())
		.def(py::init<int, int>(),
			 py::arg("nphot"), py::arg("modes"))

		.def("__str__", &basisStr)
		.def("__repr__", &basisRepr)

		.def("__len__", [](const Basis &b) { return b.size(); },
			 "Returns the number of Fock states in the basis.")

		.def("__delitem__", (void (*)(Basis &, const Fock &)) &eraseKey,
			 "Removes the Fock 'f' from the basis.",
			 py::arg("f"))

		.def("__iter__", [](Basis &b) {
				 return py::make_iterator(b.begin(), b.end());
			 },
			 "Returns corresponding iterator object.",
			 py::keep_alive<0, 1>())

		.def("__reversed__", [](Basis &b) {
				 return py::make_iterator(b.rbegin(), b.rend());
			 },
			 "Returns corresponding reverse iterator object.",
			 py::keep_alive<0, 1>())

		.def("__contains__", [](const Basis &b, const Fock &f) {
				 return b.find(f) != b.end();
			 },
			 "Tests whether the Fock 'f' is in the basis.",
			 py::arg("f"))

		.def("add", [](Basis &b, const Fock &f) {
				 b.insert(f);
			 },
			 "Adds the Fock state 'f' to the basis.",
			 py::arg("f"))

		.def("remove", (void (*)(Basis &, const Fock &)) &eraseKey,
			 "Removes the Fock state 'f' from the basis. Throws KeyError if 'f'\n"
			 "does not exist",
			 py::arg("f"))

		.def("discard", [](Basis &b, const Fock &f) {
				 b.erase(f);
			 },
			 "Removes the Fock state 'f' from the basis if it is present.",
			 py::arg("f"))

		.def("clear", [](Basis &b) { b.clear(); },
			 "Remove all elements from the basis.")

		.def(py::self + py::self,
			 "Returns a basis which is a union of Fock states from both bases.")
		.def(py::self += py::self,
			 "Effectively equivalent to self = self + b.",
			 py::arg("b"))

		.def(py::self * py::self,
			 "Calculates a tensor product of two bases.\n"
			 "Returns a basis consisting of all possible elementwise tensor\n"
			 "products of elements of the bases.")
		.def(py::self *= py::self,
			 "Effectively equivalent to self = self * b.",
			 py::arg("b"))

		.def("generate_basis", &Basis::generateBasis,
			 "Generates a basis of all possible Fock states with 'modes' modes\n"
			 "and containing 'nphot' photons.",
			 py::arg("nphot"), py::arg("modes"), py::arg("head") = Fock())

		.def("postselect", &Basis::postselect,
			 "Returns a postselected basis after observing ancilla 'anc'.\n"
			 "Ancilla is assumed to occupy the first modes.",
			 py::arg("anc"))

		.def("apply_function", [](const Basis &b, const FockAmpFunction &f, pyexecution exec_policy) {
				 State s;
				 switch(exec_policy)
				 {
				 case pyexecution::seq:
					 s =  b.applyFunction<execution::Seq>(f);
					 break;
				 case pyexecution::par:
					 s =  b.applyFunction<execution::Par>(f);
					 break;
				 }
				 return s;
			 },
			 "Applies a function 'func' to all Fock states of the basis to\n"
			 "compute a corresponding amplitude and returns the corresponding\n"
			 "state. The 'func' should take a fock object as an argument and\n"
			 "return a complex number representing its amplitude.",
			 py::arg("func"), py::arg("exec_policy") = pyexecution::seq)

		.def("as_set", [](const Basis &b) {
				 py::set s;
				 for(const auto &f: b)
					 s.add(f);
				 return s;
			 },
			 "Returns Python set object that represents the basis.")
	;

	// State
	py::class_<State>(m, "State")
		.def(py::init<>())
		.def(py::init<const State &>())
		.def(py::init<const Fock &>())
		.def(py::init<const State::Map &>())
		.def(py::init<const Basis &>())

		.def("__str__", &stateStr)
		.def("__repr__", &stateRepr)

		.def("__len__", [](const State &s) { return s.size(); },
			 "Returns number of Focks with a specified amplitude in the state.")

		.def("__getitem__", [](const State &s, const Fock &f) {
				 auto iter = s.find(f);
				 if(iter == s.end())
					 throw py::key_error();
				 return iter->second;
			 },
			 "Returns an amplitude corresponding to the Fock 'f'.",
			 py::arg("f"))

		.def("__missing__", [](State &s) { return State::Value(0.); })

		.def("__setitem__", [](State &s, const Fock &f, State::Value amp) {
				 if(amp == 0.)
					 s.erase(f);
				 else
					 s[f] = amp;
				 return amp;
			 },
			 "Sets amplitude of the Fock 'f' to 'amp'. If 'amp' = 0, then\n"
			 "deletes the Fock 'f' from the state.",
			 py::arg("f"), py::arg("amp"))

		.def("__delitem__", (void (*)(State &, const Fock &)) &eraseKey,
			 "Removes the Fock 'f' from the state. Equivalent to self[f] = 0.",
			 py::arg("f"))

		.def("__iter__", [](State &s) {
				 return py::make_key_iterator(s.begin(), s.end());
			 },
			 "Returns corresponding iterator object.",
			 py::keep_alive<0, 1>())

		.def("__reversed__", [](State &s) {
				 return py::make_key_iterator(s.rbegin(), s.rend());
			 },
			 "Returns corresponding reverse iterator object.",
			 py::keep_alive<0, 1>())

		.def("__contains__", [](const State &s, const Fock &f) {
				 return s.find(f) != s.end();
			 },
			 "Tests whether the state has a specified amplitude corresponding\n"
			 "to the Fock 'f'.",
			 py::arg("f"))

		.def("clear", [](State &s) { s.clear(); },
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

		.def(py::self * Complex(),
			 "Multiplies a state by a complex number.")
		.def(py::self *= Complex(),
			 "Effectively equivalent to self = self * z.",
			 py::arg("z"))

		.def(py::self / Complex(),
			 "Divides a state by a complex number.")
		.def(py::self /= Complex(),
			 "Effectively equivalent to self = self / z.",
			 py::arg("z"))

		.def("norm", &State::norm,
			 "Returns norm of the state.")

		.def("normalize", &State::normalize,
			 "Normalizes the state to have unit norm.")

		.def("dot", &State::dot,
			 "Calculates a dot (scalar) product.")

		.def("postselect", (State (State::*)(const Fock&) const) &State::postselect,
			 "Returns postselected state for ancilla 'anc'. The ancilla is\n"
			 "assumed to occupy the first modes.",
			 py::arg("anc"))

		.def("postselect", (std::map<Fock, State> (State::*)(int) const) &State::postselect,
			 "Calculates postselected states for all possible ancillas that\n"
			 "occupy the first 'nmodes' modes. Returns a dictionary in the\n"
			 "following format:\n"
			 "{anc1: postselected_state1, anc2: postselected_state2, ...},\n"
			 "where anc1, anc2, ... have 'nmodes' modes.\n"
			 "This function is useful when postselection for all possible\n"
			 "ancillas is required. It is faster than calling postselect(anc)\n"
			 "in cycle for different 'anc'.",
			 py::arg("nmodes"))

		.def("postselect", (std::map<Fock, State> (State::*)(const Basis &) const) &State::postselect,
			 "Calculates postselected states for all ancillas specified by\n"
			 "'basis'. Returns a dictionary in the following format:\n"
			 "{anc1: postselected_state1, anc2: postselected_state2, ...}.\n"
			 "This function is useful when postselection for an array of\n"
			 "ancillas is required. It is faster than calling postselect(anc)\n"
			 "in cycle for different 'anc'.",
			 py::arg("basis"))

		.def("get_basis", &State::getBasis,
			 "Returns basis of the state.")

		.def("set_basis", &State::setBasis,
			 "Sets basis of the state.",
			 py::arg("basis"))

		.def("get_amplitudes", &State::getAmplitudes,
			 "Returns an array of state amplitudes.")

		.def("set_amplitudes", [](State &s, const std::vector<Complex> &amps, pyexecution exec_policy) {
				switch(exec_policy)
				{
				case pyexecution::seq:
					s.setAmplitudes<execution::Seq>(amps);
					break;
				case pyexecution::par:
					s.setAmplitudes<execution::Par>(amps);
					break;
				}
			 },
			 "Sets amplitudes of the state according to the array 'amps'.",
			 py::arg("amps"), py::arg("exec_policy") = pyexecution::seq)

		.def("set_amplitudes", [](State &s, const FockAmpFunction &f, pyexecution exec_policy) {
				 switch(exec_policy)
				 {
				 case pyexecution::seq:
					 s.setAmplitudes<execution::Seq>(f);
					 break;
				 case pyexecution::par:
					 s.setAmplitudes<execution::Par>(f);
					 break;
				 }
			  },
			  "Sets amplitudes of the state using a function 'func'.",
			  py::arg("func"), py::arg("exec_policy") = pyexecution::seq)

		.def("as_dict", [](const State &s) {
				 py::dict d;
				 for(const auto &elem: s)
					d[py::make_tuple(elem.first)] = elem.second;
				 return d;
			 },
			 "Returns Python dict object that represents the state.")
	;

	m
	.def("__mul__", (State (*)(Complex, const State &)) &operator*,
		 "Multiplies a state by a complex number.",
		 py::arg("z"), py::arg("s"))
	.def("dot", &dot,
		 "Calculates a dot (scalar) product of two states.",
		 py::arg("s1"), py::arg("s2"))
	;

	// Implicit conversions
	py::implicitly_convertible<Fock::Vector, Fock>();
	py::implicitly_convertible<Basis::Set, Basis>();
	py::implicitly_convertible<State::Map, State>();
	py::implicitly_convertible<Fock, State>();

	// Circuit
	py::class_<Circuit>(m, "Circuit")
		.def(py::init<>())

		.def_property("input_state",
			 &Circuit::getInputState,
			 &Circuit::setInputState,
			 "Input state.")

		.def_property("output_basis",
			 &Circuit::getOutputBasis,
			 &Circuit::setOutputBasis,
			 "A basis of Fock states used to compute the output state. It may\n"
			 "either contain the full Fock state basis of the system or be\n"
			 "composed only of the combination of the possible states in the\n"
			 "logical and ancilla subsystems. Such a combination eliminates\n"
			 "the necessity to compute the probabilities for the output states\n"
			 "which won't be present in the system. For example, if the total\n"
			 "number of photons in the system is 4, the ancilla subsystem is\n"
			 "set to have 2 photons, and the logical subsystem is also set to\n"
			 "have 2 photons, the states, where 3 or 4 photons populate either\n"
			 "ancilla or logical subsystem are irrelevant to the problem being\n"
			 "solved.")

		.def_property("unitary",
			 &Circuit::getUnitary,
			 &Circuit::setUnitary,
			 "Unitary matrix representing a transformation of creation\n"
			 "operators of photons in modes.")

		.def("output_state", [](Circuit &c, pyexecution exec_policy) {
				 State s;
				 switch(exec_policy)
				 {
				 case pyexecution::seq:
					 s = c.outputState<execution::Seq>();
					 break;
				 case pyexecution::par:
					 s = c.outputState<execution::Par>();
					 break;
				 }
				 return s;
			 },
			 "Calculates output state.",
			 py::arg("exec_policy") = pyexecution::seq)
	;

	// Circuit design
	m
	.def("clements_design", [](Matrix Min, const Point &x, const Point &y) {
			 clementsDesign(Min, x, y);
			 return Min;
		 },
		 "Returns a unitary NxN matrix 'M' according to the Clements design.\n"
		 "Min - a diagonal unitary matrix of size NxN.\n"
		 "x - an array of N*(N-1) phase-shift parameters, such that even\n"
		 "elements are phase shifts before the beam splitters, and odd ones\n"
		 "are phase shifts between the beam splitters. All pairs of parameters\n"
		 "go in reverse column-wise enumeration order.\n"
		 "y - an array of N*(N-1) beam splitters angle defects, such that\n"
		 "even elements are defects of the first splitters, and odd ones are\n"
		 "defects of the second splitters. All pairs of parameters go in\n"
		 "reverse column-wise enumeration order.",
		 py::arg("Min"), py::arg("x"), py::arg("y"))

	.def("clements_design", (Matrix (*)(const Point &, const Point &)) &clementsDesign,
		 "Effectively equivalent to clements_design(Min, x, y) with Min being\n"
		 "the identity matrix.",
		 py::arg("x"), py::arg("y"))

	.def("clements_design", [](Matrix Min, const Point &x) {
			 clementsDesign(Min, x);
			 return Min;
		 },
		 "Effectively equivalent to clements_design(Min, x, y) with y = 0.\n"
		 "This function is provided for optimization.",
		 py::arg("Min"), py::arg("x"))

	.def("clements_design", (Matrix (*)(const Point &)) &clementsDesign,
		 "Effectively equivalent to clements_design(Min, x) with Min being\n"
		 "the identity matrix.",
		 py::arg("x"))

	.def("get_clements_design", [](Matrix M, Real eps) {
			 Point x;
			 getClementsDesign(M, x, eps);
			 return std::pair<Point, Matrix>(x, M);
		 },
		 "Calculates phase-shift coefficients 'x' for a unitary matrix 'M'\n"
		 "according to the Clements design. This function is inverse of\n"
		 "clements_design(x). Returns a pair ('x', 'Mout'), where 'Mout' is a\n"
		 "diagonal unitary matrix.\n"
		 "'eps' is a precision for a unitarity test of the input matrix 'M'.\n"
		 "If 'eps' is negative then no tests are performed.",
		 py::arg("M"), py::arg("eps") = defaultEpsilon)
	;

#if PYBIND11_VERSION_MAJOR <= 2 && PYBIND11_VERSION_MINOR < 2
	return m.ptr();
#endif
}

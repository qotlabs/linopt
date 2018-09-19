#include <iterator>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
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

// fock to list conversion functions
//py::list fock_to_list(const fock &f)
//{
//	py::list l;
//	for (unsigned int i = 0; i < f.size(); ++i)
//		l.append(f[i]);
//	return l;
//}

// state to dict conversion
//py::dict state_to_dict(const state &s)
//{
//	py::dict d;
//	for(auto iter = s.begin(); iter != s.end(); ++iter)
//	{
//		d[py::make_tuple(fock_to_list(iter -> first))] = iter -> second;
//		//d.second = iter -> second;
//		//d[iter -> first] = iter->second;
//	}
//	return d;
//}

// basis to list conversion
//py::list basis_to_list(const basis &b)
//{
//	py::list s;
//	for(auto iter: b)
//		s.append(iter);
//	return s;
//}

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

		.def(py::self < py::self)
		.def(py::self <= py::self)
		.def(py::self == py::self)
		.def(py::self != py::self)
		.def(py::self >= py::self)
		.def(py::self > py::self)

		.def("__len__", [](const fock &f) { return f.size(); })

		.def("__iter__", [](fock &f) { return py::make_iterator(f.begin(), f.end()); },
			 py::keep_alive<0, 1>())

		.def("__getitem__", [](const fock &f, int i) { return f[i]; })

		.def("__setitem__", [](fock &f, int i, int val) { return f[i] = val; })

		.def("total", &fock::total,
			 "Calculates the total number of photons in all modes.")

		.def("prod_fact", &fock::prod_fact,
			 "Calculates product of factorials of ocupation numbers.")

		.def(py::self * py::self,
			 "Calculates a tensor product of two Fock states.")
		.def(py::self *= py::self,
			 "Effectively equivalent to f1 = f1 * f2.")

//		.def_property("as_list", &fock_to_list, &list_to_fock)
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

		.def("__len__", [](const basis &b) { return b.size(); })

		.def("__iter__", [](basis &b) { return py::make_iterator(b.begin(), b.end()); },
			 py::keep_alive<0, 1>())

		.def(py::self + py::self,
			 "Returns a basis which is a union of Fock states from both bases.")
		.def(py::self += py::self,
			 "Effectively equivalent to b1 = b1 + b2.")

		.def(py::self * py::self,
			 "Calculates a tensor product of two bases.\n"
			 "Returns a basis consisting of all possible elementwise tensor "
			 "products of elements of the bases.")
		.def(py::self *= py::self,
			 "Effectively equivalent to b1 = b1*b2.")

		.def("postselect", &basis::postselect,
			 "Returns a postselected basis after observing ancilla 'anc'.\n"
			 "Ancilla is assumed to occupy the first modes.",
			 py::arg("anc"))

		.def("apply_func", &basis::apply_func,
			 "Applies a function 'func' to all Fock states of the basis to "
			 "compute corresponding amplitude and returns the corresponding "
			 "state.",
			 py::arg("func"))

//		.def_property("as_list", &basis_to_list, &list_to_basis)
	;

	// State
	py::class_<state>(m, "state")
		.def(py::init<>())
		.def(py::init<const state &>())
		.def(py::init<const state::map_class &>())

		.def("__str__", &state_str)
		.def("__repr__", &state_repr)

		.def(py::self + py::self,
			 "Adds two states, i.e., calculates their superposition.")
		.def(py::self += py::self,
			 "Effectively equivalent to s1 = s1 + s2.")

		.def(py::self - py::self,
			 "Substracts two states.")
		.def(py::self -= py::self,
			 "Effectively equivalent to s1 = s1 - s2.")

		.def(-py::self,
			 "Negates amplitudes of the state.")

		.def(py::self * py::self,
			 "Returns a tensor product of two states.")
		.def(py::self *= py::self,
			 "Effectively equivalent to s1 = s1 * s2.")

		.def(py::self * complex_type(),
			 "Multiplies a state by a complex number.")
		.def(py::self *= complex_type(),
			 "Effectively equivalent to s = s * z.")

		 .def(py::self / complex_type(),
			  "Divides a state by a complex number.")
		 .def(py::self /= complex_type(),
			  "Effectively equivalent to s = s / z.")

		.def("norm", &state::norm,
			 "Returns norm of the state.")

		.def("normalize", &state::normalize,
			 "Normalizes the state to have unit norm.")

		.def("dot", &state::dot,
			 "Calculates a dot (scalar) product.")

		.def("postselect", (state (state::*)(const fock&) const) &state::postselect,
			 "",
			 py::arg("anc"))

		.def("postselect", (std::map<fock, state> (state::*)(int) const) &state::postselect,
			 "",
			 py::arg("modes"))

		.def("get_basis", &state::get_basis,
			 "Returns basis of the state.")

//		.def_property("as_dict", &state_to_dict, &dict_to_state)
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

		.def("output_state", &circuit::output_state)

		.def_property("input_state",
					  py::cpp_function((const fock& (circuit::*)() const) &circuit::input_state, py::return_value_policy::copy),
					  &circuit::set_input_state)

		.def_property("output_basis",
					  py::cpp_function((const basis& (circuit::*)() const) &circuit::output_basis, py::return_value_policy::copy),
					  &circuit::set_output_basis)

		.def_property("unitary",
					  py::cpp_function((const matrix_type& (circuit::*)() const) &circuit::unitary, py::return_value_policy::copy),
					  &circuit::set_unitary)
	;

#if PYBIND11_VERSION_MAJOR <= 2 && PYBIND11_VERSION_MINOR < 2
	return m.ptr();
#endif
}

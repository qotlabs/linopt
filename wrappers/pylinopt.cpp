#include <iterator>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
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

//void list_to_fock(fock &f, const py::list &l)
//{
//	f = {};
//	if(l.size() > 0)
//	{
//		for (unsigned int i = 0; i < l.size(); ++i)
//			f.push_back(py::cast<double>(l[i]));
//	}
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

//void dict_to_state(state &s, const py::dict &d)
//{
//	//py::list keys = d.keys();
//	for(auto item : d)
//	{
//		s[py::cast<fock>(item.first)] = py::cast<complex_type>(item.second);
//	}
//}

// basis to list conversion
//py::list basis_to_list(const basis &b)
//{
//	py::list s;
//	for(auto iter: b)
//		s.append(iter);
//	return s;
//}

//void list_to_basis(basis &b, const py::list &l)
//{
//	b = {};
//	if(l.size() > 0)
//	{
//		for (unsigned int i = 0; i < l.size(); ++i)
//			b.insert(py::cast<fock>(l[i]));
//	}
//}

//matrix_type exp_hermite(const py::array_t<double> array)
//{
//	point vec(array.size());
//	std::memcpy(vec.data(),array.data(),array.size()*sizeof(double));
//	return exp_hermite_parametrization(vec);
//}

//matrix_type hurwitz(const py::array_t<double> array)
//{
//	point vec(array.size());
//	std::memcpy(vec.data(),array.data(),array.size()*sizeof(double));
//	return hurwitz_parametrization(vec);
//}

//void exp_hermite(matrix_type &M, const py::array_t<double> array)
//{
//	M = exp_hermite(array);
//}

//void hurwitz(matrix_type &M, const py::array_t<double> array)
//{
//	M = hurwitz(array);
//}

//basis& generate_basis(basis &b, const int nphot, const int modes)
//{
//	return b.generate_basis(nphot, modes);
//}

//PYBIND11_MAKE_OPAQUE(fock::vector)
//PYBIND11_MAKE_OPAQUE(point)

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

//	py::bind_vector<fock::vector>(m, "fock_vector");
//	py::bind_vector<point>(m, "point");

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
		.def("total", &fock::total)
		.def("prod_fact", &fock::prod_fact)
//		.def("__mul__", &fock::operator*)
//		.def("__imul__", &fock::operator*=, py::return_value_policy::copy)
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
//		.def("__add__", &basis::operator+)
//		.def("__iadd__", &basis::operator+=, py::return_value_policy::copy)
//		.def("__mul__", &basis::operator*)
//		.def("__imul__", &basis::operator*=, py::return_value_policy::copy)
//		.def("generate_basis", &generate_basis, py::return_value_policy::copy)
//		.def("postselect", &basis::postselect)
//		.def("apply_func", &basis::apply_func)
//		.def_property("as_list", &basis_to_list, &list_to_basis)
	;

	// State
	py::class_<state>(m, "state")
		.def(py::init<>())
		.def(py::init<const state &>())
		.def(py::init<const state::map_class &>())
		.def("__str__", &state_str)
		.def("__repr__", &state_repr)
//		.def("__add__", (state (state::*)(const state &) const) &state::operator+)
//		.def("__iadd__", (state& (state::*)(const state &)) &state::operator+=, py::return_value_policy::copy)
//		.def("__sub__", (state (state::*)() const) &state::operator-)
//		.def("__sub__", (state (state::*)(const state &) const) &state::operator-)
//		.def("__isub__", (state& (state::*)(const state &)) &state::operator-=, py::return_value_policy::copy)
//		.def("__mul__", (state (state::*)(const state &) const) &state::operator*)
//		.def("__mul__", (state (state::*)(complex_type) const) 	&state::operator*)
//		.def("__imul__", (state& (state::*)(const state &)) &state::operator*=, py::return_value_policy::copy)
//		.def("__imul__", (state& (state::*)(complex_type)) &state::operator*=, py::return_value_policy::copy)
//		.def("__div__", (state (state::*)(complex_type) const) &state::operator/)
//		.def("__itruediv__", (state& (state::*)(complex_type)) &state::operator/=, py::return_value_policy::copy)
//		.def("norm", &state::norm)
//		.def("normalize", &state::normalize, py::return_value_policy::copy)
//		.def("dot", &state::dot)
//		.def("postselect", (state (state::*)(const fock&) const) &state::postselect)
//		.def("get_basis", &state::get_basis)
//		.def_property("as_dict", &state_to_dict, &dict_to_state)
	;

	// Circuit
//	py::class_<circuit>(m, "circuit")
//		.def(py::init<>())
//		.def("output_state", &circuit::output_state)
//		.def_property("input_state",
//					  py::cpp_function((const fock& (circuit::*)() const) &circuit::input_state, py::return_value_policy::copy),
//					  &circuit::set_input_state)
//		.def_property("output_basis",
//					  py::cpp_function((const basis& (circuit::*)() const) &circuit::output_basis, py::return_value_policy::copy),
//					  &circuit::set_output_basis)
//		.def_property("unitary",
//					  py::cpp_function((const matrix_type& (circuit::*)() const) &circuit::unitary, py::return_value_policy::copy),
//					  &circuit::set_unitary)
//	;

	py::implicitly_convertible<fock::vector_class, fock>();

#if PYBIND11_VERSION_MAJOR <= 2 && PYBIND11_VERSION_MINOR < 2
	return m.ptr();
#endif
}

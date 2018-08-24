#include <iterator>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include "../lib/linopt.h"

using namespace linopt;
namespace py = pybind11;

// __repr__ and __str__ methods for class "fock"
std::string fock_str(fock& _state)
{   
	std::stringstream buffer;

	if(_state.size() > 0)
	{
	   buffer << "[";
	   std::copy(_state.begin(), _state.end()-1, std::ostream_iterator<int>(buffer, ", "));
	   buffer << _state.back() << "]"; 
	}
	else
	{
		buffer << "Empty fock";
	}
	return buffer.str();
}

std::string fock_repr(fock& _state)
{   
	std::stringstream buffer;

	if(_state.size() > 0)
	{
		buffer << "Fock state: [";
		std::copy(_state.begin(), _state.end()-1, std::ostream_iterator<int>(buffer, ", "));
		buffer << _state.back() << "]";
	}
	else
	{
		buffer << "Empty fock";
	}
	
	return buffer.str();
}

// __repr__ and __str__ methods for class "basis"
std::string basis_str(basis& bas)
{   
	std::stringstream buffer;

	if(bas.size() > 0)
	{
		buffer << "{ ";
		for (const fock& el : bas)
		{
			buffer << el << ' ';
		}
		buffer << " }"; 
	}
	else
	{
		buffer << "Empty basis";
	}
	
	return buffer.str();
}

std::string basis_repr(basis& bas)
{   
	std::stringstream buffer;

	if(bas.size() > 0)
	{
		buffer << "Fock state basis: { ";
		for (const fock& el : bas)
		{
			buffer << el << ' ';
		}
		buffer << " }";   
	}
	else
	{
		buffer << "Empty basis";
	}
	
	return buffer.str();
}

// __repr__ and __str__ methods for class "state"
std::string state_str(state& sta)
{
	std::stringstream buffer;

	if(sta.size() > 0)
	{
		buffer << "{ ";
		for(auto it = sta.cbegin(); it != sta.cend(); ++it)
		{
			buffer << it->first << ": " << it->second;
		}
		buffer << " }";    
	}
	else
	{
		buffer << "Empty state";
	}
	
	return buffer.str();
}

std::string state_repr(state& sta)
{
	std::stringstream buffer;

	if(sta.size() > 0)
	{
		buffer << "Fock state coefficients: { ";
		for(auto it = sta.cbegin(); it != sta.cend(); ++it)
		{
			buffer << it->first << ": " << it->second << " ";
		}
		buffer << " }";    
	}
	else
	{
		buffer << "Empty state";
	}
	
	return buffer.str();
}

// fock to list conversion functions
py::list fock_to_list(const fock &f)
{
	py::list l;
	for (unsigned int i = 0; i < f.size(); ++i)
		l.append(f[i]);
	return l;
}

void list_to_fock(fock &f, const py::list &l)
{
	f = {};
	if(l.size() > 0)
	{
		for (unsigned int i = 0; i < l.size(); ++i)
			f.push_back(py::cast<double>(l[i]));
	}
}

// state to dict conversion
py::dict state_to_dict(const state &s)
{
	py::dict d;
	for(auto iter = s.begin(); iter != s.end(); ++iter)
	{
		d[py::make_tuple(fock_to_list(iter -> first))] = iter -> second;
		//d.second = iter -> second;
		//d[iter -> first] = iter->second;
	}
	return d;
}

void dict_to_state(state &s, const py::dict &d)
{
	//py::list keys = d.keys();
	for(auto item : d)
	{
		s[py::cast<fock>(item.first)] = py::cast<complex_type>(item.second);
	}
}

// basis to list conversion
py::list basis_to_list(const basis &b)
{
	py::list s;
	for(auto iter : b)
		s.append(iter);
	return s;
}

void list_to_basis(basis &b, const py::list &l)
{
	b = {};
	if(l.size() > 0)
	{
		for (unsigned int i = 0; i < l.size(); ++i)
			b.insert(py::cast<fock>(l[i]));
	}
}

matrix_type exp_hermite(const py::array_t<double> array)
{
	point vec(array.size());
	std::memcpy(vec.data(),array.data(),array.size()*sizeof(double));
	return exp_hermite_parametrization(vec);
}

matrix_type hurwitz(const py::array_t<double> array)
{
	point vec(array.size());
	std::memcpy(vec.data(),array.data(),array.size()*sizeof(double));
	return hurwitz_parametrization(vec);
}

void exp_hermite(matrix_type &M, const py::array_t<double> array)
{
	M = exp_hermite(array);
}

void hurwitz(matrix_type &M, const py::array_t<double> array)
{
	M = hurwitz(array);
}
// "thin wrappers" for methods with default arguments from class unitary_matrix
// bool is_column_unitary_noargs(const matrix_type &M)
// { 
// 	return u.is_column_unitary();
// }
// bool is_row_unitary_noargs(const matrix_type &M)
// { 
// 	return u.is_row_unitary();
// }
// bool is_unitary_noargs(const matrix_type &M)
// { 
// 	return u.is_unitary();
// }

basis& generate_basis(basis &b, const int nphot, const int modes)
{
	return b.generate_basis(nphot, modes);
}

PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<complex_type>);

PYBIND11_PLUGIN(pylinopt)
{   
	// Py_Initialize();
	// np::initialize();
	// Register interable conversions.
	// iterable_converter()
	// Built-in type.
	// .from_python<std::vector<double> >()
	// ;

	// class_<vector_type>("vector_type")
	// 	.def(init<>())
	// 	.def(VectorVisitor<vector_type>())
	// ;

	// class_<matrix_type>("matrix_type")
	// 	.def(init<>())
	// 	.def(MatrixVisitor<matrix_type>())
	// ;

	// def("float2str", &doubleToShortest, (arg("f"), arg("pad")=0), "Return the shortest string representation of *f* which will is equal to *f* when converted back to float. This function is only useful in Python prior to 3.0; starting from that version, standard string conversion does just that.");

	py::module m("pylinopt", "documentation string");

	// exposing matrix functions
	m.def("permanent", &permanent);
	m.def("hurwitz", (matrix_type (*)(const py::array_t<double> array)) &hurwitz);
	m.def("hurwitz", (void (*)(matrix_type &M, const py::array_t<double> array)) &hurwitz);
	m.def("exp_hermite", (matrix_type (*)(const py::array_t<double> array)) &exp_hermite);
	m.def("exp_hermite", (void (*)(matrix_type &M, const py::array_t<double> array)) &exp_hermite);

	// std::vector<int> class declaration (class fock inherits from this class)
	py::bind_vector<std::vector<int>>(m, "VectorInt");

	// expose class fock
	py::class_< fock, std::vector<int> >(m, "fock")
		.def(py::init<>())
		.def("__str__", &fock_str)
		.def("__repr__", &fock_repr)
		.def("total", &fock::total)
		.def("prod_fact", &fock::prod_fact)
		.def("__mul__", &fock::operator*)
		.def("__imul__", &fock::operator*=, py::return_value_policy::copy)
		.def_property("as_list", &fock_to_list, &list_to_fock)
	;

	py::bind_map<std::map<fock, complex_type>>(m, "MapFockComplex");

	// expose class state
	py::class_< state, std::map<fock, complex_type> >(m, "state")
		.def(py::init<>())
		.def(py::init<const state&>())
		.def("__str__", &state_str)
		.def("__repr__", &state_repr)
		.def("__add__", (state (state::*)(const state &) const) &state::operator+)
		.def("__iadd__", (state& (state::*)(const state &)) &state::operator+=, py::return_value_policy::copy)
		.def("__sub__", (state (state::*)() const) &state::operator-)
		.def("__sub__", (state (state::*)(const state &) const) &state::operator-)
		.def("__isub__", (state& (state::*)(const state &)) &state::operator-=, py::return_value_policy::copy)
		.def("__mul__", (state (state::*)(const state &) const) &state::operator*)
		.def("__mul__", (state (state::*)(complex_type) const) 	&state::operator*)
		.def("__imul__", (state& (state::*)(const state &)) &state::operator*=, py::return_value_policy::copy)
		.def("__imul__", (state& (state::*)(complex_type)) &state::operator*=, py::return_value_policy::copy)
		.def("__div__", (state (state::*)(complex_type) const) &state::operator/)
		.def("__itruediv__", (state& (state::*)(complex_type)) &state::operator/=, py::return_value_policy::copy)
		.def("norm", &state::norm)
		.def("normalize", &state::normalize, py::return_value_policy::copy)
		.def("dot", &state::dot)
		.def("postselect", &state::postselect)
		.def("get_basis", &state::get_basis)
		.def_property("as_dict", &state_to_dict, &dict_to_state)
	;

	py::class_< basis >(m, "basis")
		.def(py::init<>())
		.def(py::init<const basis&>())
		.def("__str__", &basis_str)
		.def("__repr__", &basis_repr)
		.def("__add__", &basis::operator+)
		.def("__iadd__", &basis::operator+=, py::return_value_policy::copy)
		.def("__mul__", &basis::operator*)
		.def("__imul__", &basis::operator*=, py::return_value_policy::copy)
		.def("generate_basis", &generate_basis, py::return_value_policy::copy)
		.def("postselect", &basis::postselect)
		.def("apply_func", &basis::apply_func)
		.def_property("as_list", &basis_to_list, &list_to_basis)
	;

	// overloaded method bindings for class circuit
	// const matrix_type& (circuit::*const_unitary)() const = &circuit::unitary;
	// const fock& (circuit::*const_input_state)() const = &circuit::input_state;
	// const basis& (circuit::*const_output_basis)() const = &circuit::output_basis;

	py::class_< circuit >(m, "circuit")
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
	return m.ptr();
}

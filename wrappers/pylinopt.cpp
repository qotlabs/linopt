#include <iterator>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// #include <boost/python.hpp>
// #include <boost/python/numpy.hpp>
// #include <boost/python/suite/indexing/vector_indexing_suite.hpp>
// #include <boost/python/suite/indexing/map_indexing_suite.hpp>
// #include <boost/python/suite/indexing/indexing_suite.hpp>
// #include "set_indexing_suite.h"

#include "../lib/states.h"
#include "../lib/matrix.h"
#include "../lib/circuit.h"
#include "../lib/cost_functor.h"

// #include "minieigen/src/common.hpp"
// #include "minieigen/src/visitors.hpp"

// using namespace boost::python;
using namespace linopt;
// namespace np = boost::python::numpy;
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
// py::list fock_to_list(const fock &f)
// {
// 	py::list l;
// 	for (unsigned int i = 0; i < f.size(); ++i)
// 		l.append(f[i]);
// 	return l;
// }

// void list_to_fock(fock &f, const py::list &l)
// {
// 	f = {};
// 	if(len(l) > 0)
// 	{
// 		for (int i = 0; i < len(l); ++i)
// 			f.push_back(py::cast<double>(l[i]));
// 	}
// }

// state to dict conversion
// py::dict state_to_dict(const state &s)
// {
// 	py::dict d;
// 	for(auto iter = s.begin(); iter != s.end(); ++iter)
// 	{
// 		d.first = iter -> first;
// 		d.second = iter -> second;
// 		//d[iter -> first] = iter->second;
// 	}
// 	return d;
// }

// void dict_to_state(state &s, const py::dict &d)
// {
// 	//py::list keys = d.keys();
// 	for(auto item : d)
// 	{
// 		s[py::cast<fock>(item.first)] = py::cast<complex_type>(item.second);
// 	}
// }

//basis to list conversion
// py::list basis_to_list(const basis &b)
// {
// 	py::list s;
// 	for(auto iter : b)
// 		s.append(iter);
// 	return s;
// }

// void list_to_basis(basis &b, const py::list &l)
// {
// 	b = {};
// 	if(len(l) > 0)
// 	{
// 		for (int i = 0; i < len(l); ++i)
// 			b.insert(py::cast<fock>(l[i]));
// 	}
// }

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
basis& generate_basis(const basis &b, const int nphot, const int modes)
{
	return b.generate_basis( nphot, modes);
}

void (*hurwitz_twoargs)(matrix_type &M, const point &x) = &hurwitz_parametrization;
matrix_type (*hurwitz_onearg)(const point &x) = &hurwitz_parametrization;
void (*exp_hermite_twoargs)(matrix_type &M, const point &x) = &exp_hermite_parametrization;
matrix_type (*exp_hermite_onearg)(const point &x) = &exp_hermite_parametrization;

// matrix_type hurwitz_parametrization_wrap(const point &x)
// {
// 	return hurwitz_parametrization(x);
// }

// matrix_type exp_hermite_parametrization_wrap(const point &x)
// {
// 	return exp_hermite_parametrization(x);
// }

// the code below is borrowed from stackoverflow answer 
// https://stackoverflow.com/questions/15842126/feeding-a-python-list-into-a-function-taking-in-a-vector-with-boost-python

/// @brief Type that allows for registration of conversions from
///        python iterable types.
// struct iterable_converter
// {
//   /// @note Registers converter from a python interable type to the
//   ///       provided type.
//   template <typename Container>
//   iterable_converter&
//   from_python()
//   {
//     converter::registry::push_back(
//       &iterable_converter::convertible,
//       &iterable_converter::construct<Container>,
//       boost::python::type_id<Container>());

//     // Support chaining.
//     return *this;
//   }

//   /// @brief Check if PyObject is iterable.
//   static void* convertible(PyObject* object)
//   {
//     return PyObject_GetIter(object) ? object : NULL;
//   }

//   /// @brief Convert iterable PyObject to C++ container type.
//   ///
//   /// Container Concept requirements:
//   ///
//   ///   * Container::value_type is CopyConstructable.
//   ///   * Container can be constructed and populated with two iterators.
//   ///     I.e. Container(begin, end)
//   template <typename Container>
//   static void construct(
//     PyObject* object,
//     converter::rvalue_from_python_stage1_data* data)
//   {
//     namespace python = boost::python;
//     // Object is a borrowed reference, so create a handle indicting it is
//     // borrowed for proper reference counting.
//     python::handle<> handle(python::borrowed(object));

//     // Obtain a handle to the memory block that the converter has allocated
//     // for the C++ type.
//     typedef python::converter::rvalue_from_python_storage<Container>
//                                                                 storage_type;
//     void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

//     typedef python::stl_input_iterator<typename Container::value_type>
//                                                                     iterator;

//     // Allocate the C++ type into the converter's memory block, and assign
//     // its handle to the converter's convertible variable.  The C++
//     // container is populated by passing the begin and end iterators of
//     // the python object to the container's constructor.
//     new (storage) Container(
//       iterator(python::object(handle)), // begin
//       iterator());                      // end
//     data->convertible = storage;
//   }
// };

// ------------------------- end of stackoverflow code ---------------------------------------

PYBIND11_MODULE(pylinopt, m)
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

	PYBIND11_MAKE_OPAQUE(std::vector<real_type>);
	PYBIND11_MAKE_OPAQUE(std::vector<complex_type>);

	// exposing matrix functions
	m.def("permanent", &permanent);
	m.def("hurwitz", &hurwitz_onearg);
	m.def("hurwitz", &hurwitz_twoargs);
	m.def("exp_hermite", &exp_hermite_onearg);
	m.def("exp_hermite", &exp_hermite_twoargs);

	// std::vector<int> class declaration (class fock inherits from this class)
	// class_<std::vector<int>>("cvector")
	// 	.def(vector_indexing_suite<std::vector<int>>())
	// ;

	// expose class fock
	py::class_< fock, bases<std::vector<int> > >(m, "fock")
		.def(init<>())
		.def("__str__", &fock_str)
		.def("__repr__", &fock_repr)
		.def("total", &fock::total)
		.def("prod_fact", &fock::prod_fact)
		.def("__mul__", &fock::operator*)
		.def("__imul__", &fock::operator*=, py::return_value_policy::copy)
		//.def_property("as_list", &fock_to_list, &list_to_fock)
	;

	// std::map<fock, complex_type> class declaration (class state inherits from this class)
	// class_< std::map<fock, complex_type> >("cmapfock")
	// 	.def(map_indexing_suite<std::map<fock, complex_type>>())
	// ;

	// overloaded operator bindings for class state
	// state (state::*unary_plus)(const state &) const 	= &state::operator+;
	// state& (state::*state_iadd)(const state &) 			= &state::operator+=;
	// state (state::*unary_minus)() const 				= &state::operator-;
	// state (state::*binary_minus)(const state &) const 	= &state::operator-;
	// state& (state::*state_isub)(const state &) 			= &state::operator-=;
	// state (state::*state_mul)(const state &) const 		= &state::operator*;
	// state (state::*complex_mul)(complex_type) const 	= &state::operator*;
	// state& (state::*state_imul)(const state &) 			= &state::operator*=;
	// state& (state::*complex_imul)(complex_type) 		= &state::operator*=;
	// state (state::*complex_div)(complex_type) const 	= &state::operator/;
	// state& (state::*complex_idiv)(complex_type) 		= &state::operator/=;

	// expose class state
	py::class_< state, bases<std::map<fock, complex_type> > >(m, "state")
		.def(init<>())
		.def(init<const state&>())
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
		//.def_property("as_dict", &state_to_dict, &dict_to_state)
	;

	// class std::set<fock> declaration (class basis inherits from this class)
	// used custom "set_indexing_suite.h" and "list_indexing_suite.h" libraries to expose std::set to python correctly
	// downloaded from: https://github.com/Microsoft/bond/blob/master/python/inc/bond/python/set_indexing_suite.h
	// class_< std::set<fock> >("csetfock")
	// 	.def(set_indexing_suite<std::set<fock>>());

	py::class_< basis, bases< std::set<fock> > >(m, "basis")
		.def(init<>())
		.def(init<const basis&>())
		.def("__str__", &basis_str)
		.def("__repr__", &basis_repr)
		.def("__add__", &basis::operator+)
		.def("__iadd__", &basis::operator+=, py::return_value_policy::copy)
		.def("__mul__", &basis::operator*)
		.def("__imul__", &basis::operator*=, py::return_value_policy::copy)
		.def("generate_basis", &generate_basis, py::return_value_policy::copy)
		.def("postselect", &basis::postselect)
		.def("apply_func", &basis::apply_func)
		//.def_property("as_list", &basis_to_list, &list_to_basis)
	;

	// overloaded method bindings for class circuit
	// const matrix_type& (circuit::*const_unitary)() const = &circuit::unitary;
	// const fock& (circuit::*const_input_state)() const = &circuit::input_state;
	// const basis& (circuit::*const_output_basis)() const = &circuit::output_basis;

	py::class_< circuit >(m, "circuit")
		.def(init<>())
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
}

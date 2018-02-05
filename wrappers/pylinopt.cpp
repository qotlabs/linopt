#include <Eigen/Dense>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/suite/indexing/indexing_suite.hpp>
#include "set_indexing_suite.h"
#include "../lib/states.h"
#include "../lib/matrix.h"
#include "../lib/chip.h"
#include "../lib/cost_functor.h"

/**** double-conversion helpers *****/
#include "../../minieigen/src/double-conversion/double-conversion.h"
#include "../../minieigen/src/visitors.hpp"

using namespace boost::python;
using namespace linopt;

// "thin wrappers" for methods with default arguments from class unitary_matrix
bool is_column_unitary_noargs(const unitary_matrix &u)
{ 
    return u.is_column_unitary();
}
bool is_row_unitary_noargs(const unitary_matrix &u)
{ 
    return u.is_row_unitary();
}
bool is_unitary_noargs(const unitary_matrix &u)
{ 
    return u.is_unitary();
}

BOOST_PYTHON_MODULE(pylinopt)
{   
    typedef double real_type;
    typedef std::complex<real_type> complex_type;
    //typedef MatrixVisitor< Eigen::Matrix< complex_type, Eigen::Dynamic, Eigen::Dynamic > > matrix_type;
    //typedef Eigen::Matrix<complex_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> matrix_type;

    class_<MatrixXcr>("MatrixX","XxX (dynamic-sized) float matrix. Constructed from list of rows (as VectorX).\n\nSupported operations (``m`` is a MatrixX, ``f`` if a float/int, ``v`` is a VectorX): ``-m``, ``m+m``, ``m+=m``, ``m-m``, ``m-=m``, ``m*f``, ``f*m``, ``m*=f``, ``m/f``, ``m/=f``, ``m*m``, ``m*=m``, ``m*v``, ``v*m``, ``m==m``, ``m!=m``.",py::init<>())
        .def(MatrixVisitor<MatrixXcr>())
    ;

    // class matrix exposing
    class_< unitary_matrix, bases<MatrixXcr> >("matrix")
        .def(init<>())
        .def("hurwitz", &unitary_matrix::hurwitz, return_value_policy<copy_non_const_reference>())
        .def("exp_hermite", &unitary_matrix::exp_hermite, return_value_policy<copy_non_const_reference>())
        .def("is_column_unitary", &unitary_matrix::is_column_unitary)
        .def("is_column_unitary", is_column_unitary_noargs)
        .def("is_row_unitary", &unitary_matrix::is_row_unitary)
        .def("is_row_unitary", is_row_unitary_noargs)
        .def("is_unitary", &unitary_matrix::is_unitary)
        .def("is_unitary", is_unitary_noargs)
        .attr("default_epsilon") = unitary_matrix::default_epsilon;

    // exposing matrix functions
    def("permanent", permanent);

    // overloaded operator bindings for class state
    state (state::*unary_minus)() const = &state::operator-;
    state (state::*binary_minus)(const state &) const = &state::operator-;
    state (state::*state_mul)(const state &) const = &state::operator*;
    state (state::*complex_mul)(complex_type) const = &state::operator*;
    state& (state::*state_imul)(const state &) = &state::operator*=;
    state& (state::*complex_imul)(complex_type) = &state::operator*=;

    // std::vector<int> class declaration (class fock inherits from this class)
    class_<std::vector<int>>("cvector")
        .def(vector_indexing_suite<std::vector<int>>());

    // fock class exposing
    class_ < fock, bases<std::vector<int> > >("fock")
		.def(init<>())
        .def("total", &fock::total)
        .def("prod_fact", &fock::prod_fact)
        .def("__mul__", &fock::operator*)
        .def("__imul__", &fock::operator*=, return_value_policy<copy_non_const_reference>());

    // std::map<fock, complex_type> class declaration (class state inherits from this class)
    class_< std::map<fock, complex_type> >("cmapfock")
        .def(map_indexing_suite<std::map<fock, complex_type>>());

    // class state exposing
    class_< state, bases<std::map<fock, complex_type> > >("state")
        .def(init<>())
        .def(init<const state&>())
        .def("__add__", &state::operator+)
        .def("__iadd__", &state::operator+=, return_value_policy<copy_non_const_reference>())
        .def("__sub__", unary_minus)
        .def("__sub__", binary_minus)
        .def("__isub__", &state::operator-=, return_value_policy<copy_non_const_reference>())
        .def("__mul__", state_mul)
        .def("__mul__", complex_mul)
        .def("__imul__", state_imul, return_value_policy<copy_non_const_reference>())
        .def("__imul__", complex_imul, return_value_policy<copy_non_const_reference>())
        .def("__div__", &state::operator/)
        .def("__idiv__", &state::operator/=, return_value_policy<copy_non_const_reference>())
        .def("norm", &state::norm)
        .def("normalize", &state::normalize, return_value_policy<copy_non_const_reference>())
        .def("dot", &state::dot)
        .def("postselect", &state::postselect)
        .def("get_basis", &state::get_basis);

    // class std::set<fock> declaration (class basis inherits from this class)
    // used custom "set_indexing_suite.h" and "list_indexing_suite.h" libraries to expose std::set to python correctly
    // downloaded from: https://github.com/Microsoft/bond/blob/master/python/inc/bond/python/set_indexing_suite.h
    // 
    class_< std::set<fock> >("csetfock")
        .def(set_indexing_suite<std::set<fock>>());

    class_< basis, bases< std::set<fock> > >("basis")
        .def(init<>())
        .def(init<const basis&>())
        .def("__add__", &basis::operator+)
        .def("__iadd__", &basis::operator+=, return_value_policy<copy_non_const_reference>())
        .def("__mul__", &basis::operator*)
        .def("__imul__", &basis::operator*=, return_value_policy<copy_non_const_reference>())
        .def("generate_basis", &basis::generate_basis, return_value_policy<copy_non_const_reference>())
        .def("postselect", &basis::postselect)
        .def("apply_func", &basis::apply_func);

    // overloaded method bindings for class chip
    unitary_matrix& (chip::*unitary)() = &chip::unitary;
    const unitary_matrix& (chip::*const_unitary)() const = &chip::unitary;
    fock& (chip::*input_state)() = &chip::input_state;
    const fock& (chip::*const_input_state)() const = &chip::input_state;
    basis& (chip::*output_basis)() = &chip::output_basis;
    const basis& (chip::*const_output_basis)() const = &chip::output_basis;

    class_< chip >("chip")
        .def(init<>())
        .def("unitary", unitary, return_value_policy<copy_non_const_reference>())
        .def("unitary", const_unitary, return_value_policy<copy_const_reference>())
        .def("input_state", input_state, return_value_policy<copy_non_const_reference>())
        .def("input_state", const_input_state, return_value_policy<copy_const_reference>())
        .def("output_basis", output_basis, return_value_policy<copy_non_const_reference>())
        .def("output_basis", const_output_basis, return_value_policy<copy_const_reference>())
        .def("output_state", &chip::output_state);

    class_< cost_functor >("cost_functor", no_init)
         .def(init<const basis&, const basis&, const fock&, const std::vector<state> >());

    class_< stanisic_functor, bases<cost_functor> >("stanisic_functor", no_init)
         .def("__call__", &stanisic_functor::operator());

    class_< log_functor, bases<cost_functor> >("log_functor", no_init)
         .def("__call__", &log_functor::operator());
}

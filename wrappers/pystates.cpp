#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/suite/indexing/indexing_suite.hpp>
#include "set_indexing_suite.h"
#include "../lib/states.h"

using namespace boost::python;
using namespace linopt;

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(sub_overloads, state::operator-, 0, 1);

BOOST_PYTHON_MODULE(pystates)
{   
    state (state::*unary_minus)() const = &state::operator-;
    state (state::*binary_minus)(const state &) const = &state::operator-;

    class_<std::vector<int>>("cvector")
        .def(vector_indexing_suite<std::vector<int>>());

    class_ < fock, bases<std::vector<int> > >("fock")
		.def(init<>())
        .def("total", &fock::total)
        .def("prod_fact", &fock::prod_fact)
        .def("__mul__", &fock::operator*)
        .def("__imul__", &fock::operator*=, return_value_policy<copy_non_const_reference>());

    class_< std::map<fock, complex_type> >("cmapfock")
        .def(map_indexing_suite<std::map<fock, complex_type>>());

    class_< state, bases<std::map<fock, complex_type> > >("state")
        .def(init<>())
        .def(init<const state&>())
        .def("__add__", &state::operator+)
        .def("__iadd__", &state::operator+=, return_value_policy<copy_non_const_reference>())

        //.def(self-self);
        //.def("__sub__", &state::operator-, sub_overloads());
        .def("__sub__", unary_minus)
        .def("__sub__", binary_minus);

        // .def("__isub__", &state::operator-=, return_value_policy<copy_non_const_reference>())
        
        // .def("__mul__", &state::operator*, args("f"))
        // .def("__mul__", &state::operator*=, args("f"), return_value_policy<copy_non_const_reference>())

    class_< std::set<fock> >("csetfock")
        .def(set_indexing_suite<std::set<fock>>());
}

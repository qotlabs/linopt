#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "../lib/states.h"

using namespace boost::python;
using namespace linopt;

BOOST_PYTHON_MODULE(pystates)
{   
    class_<std::vector<int>>("cvector")
        .def(vector_indexing_suite<std::vector<int>>());

    class_ < fock, bases<std::vector<int> > >("fock")
		.def(init<>())
        .def("total", &fock::total)
        .def("prod_fact", &fock::prod_fact)
        .def("__mul__", &fock::operator*, args( "f" ))
        .def("__imul__", &fock::operator*=, args( "f" ), return_value_policy<copy_non_const_reference>());
}

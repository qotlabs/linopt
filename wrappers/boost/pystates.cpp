#include <boost/python.hpp>
#include "../../lib/states.h"

using namespace boost::python;
using namespace linopt;

BOOST_PYTHON_MODULE(pystates)
{
    class_ < fock, bases<std::vector<int> > >("fock", no_init)
        .def("total", &fock::total)
        .def("prod_fact", &fock::prod_fact)
        .def("__mul__", &fock::operator*, args( "f" ))
        .def("__imul__", &fock::operator*=, args( "f" ), return_value_policy<copy_non_const_reference>());
}
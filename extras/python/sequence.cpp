#include "sequence.h"
#include <boost/python.hpp>

using namespace boost::python;
using namespace btl;

BOOST_PYTHON_MODULE(libbtl_py)
{
    class_<Sequence>("Sequence")
        .def("reverseComplement", &Sequence::reverseComplement)
    ;
}
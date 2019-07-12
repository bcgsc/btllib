#include "basic_sequence.h"
#include "sequence.h"
#include <boost/python.hpp>

using namespace boost::python;
using namespace btl;

BOOST_PYTHON_MODULE(btl_py)
{
    class_<BasicSequence>("BasicSequence")
        .def(self_ns::str(self)); 
    ;

    class_<Sequence, bases<BasicSequence>>("Sequence")
        .def(init<>())
        .def(init<const BasicSequence&>())
        .def(init<const Sequence&>())
        .def(init<std::string>())
        .def(init<const Sequence&, size_t, size_t>())
        .def(init<const std::string&, size_t, size_t>())
        .def(init<const char*>())
        .def(init<const char*, size_t>())
        .def(init<size_t, char>())
        .def(init<std::initializer_list<char>>())
        .def("reverseComplement", &Sequence::reverseComplement)
        .def("getReverseComplement", &Sequence::getReverseComplement)
        .def(~self)
    ;
}
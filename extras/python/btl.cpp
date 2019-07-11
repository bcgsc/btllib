#include "basic_nt_sequence.h"
#include "nt_sequence.h"
#include <boost/python.hpp>

using namespace boost::python;
using namespace btl;

BOOST_PYTHON_MODULE(btl_py)
{
    class_<BasicNtSequence>("BasicNtSequence")
        .def(self_ns::str(self)); 
    ;

    class_<NtSequence, bases<BasicNtSequence>>("NtSequence")
        .def(init<>())
        .def(init<const BasicNtSequence&>())
        .def(init<const NtSequence&>())
        .def(init<std::string>())
        .def(init<const NtSequence&, size_t, size_t>())
        .def(init<const std::string&, size_t, size_t>())
        .def(init<const char*>())
        .def(init<const char*, size_t>())
        .def(init<size_t, char>())
        .def(init<std::initializer_list<char>>())
        .def("reverseComplement", &NtSequence::reverseComplement)
        .def("getReverseComplement", &NtSequence::getReverseComplement)
        .def(~self)
    ;
}
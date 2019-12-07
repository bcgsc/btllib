%module btl

%{
#define SWIG_FILE_WITH_INIT

#include "btl/graph.hpp"
#include "btl/seq_reader.hpp"
#include "btl/seq.hpp"
#include "btl/status.hpp"
%}

%include <pyprimtypes.swg>
%include <pyopers.swg>
%include <std_common.i>
%include <cstring.i>
%include <exception.i>
%include <std_string.i>
%include <std_iostream.i>

%include "extra.i"

%include "btl/graph.hpp"
%include "btl/seq_reader.hpp"
%include "btl/seq.hpp"
%include "btl/status.hpp"

%module btllib

%{
#define SWIG_FILE_WITH_INIT

#include "btllib/seq_reader.hpp"
#include "btllib/graph.hpp"
#include "btllib/util.hpp"
#include "btllib/data_saveload.hpp"
#include "btllib/status.hpp"
#include "btllib/seq_writer.hpp"
#include "btllib/seq.hpp"
%}

%include <pyprimtypes.swg>
%include <pyopers.swg>
%include <std_common.i>
%include <cstring.i>
%include <exception.i>
%include <std_string.i>
%include <std_iostream.i>

%include "extra.i"

%include "btllib/seq_reader.hpp"
%include "btllib/graph.hpp"
%include "btllib/util.hpp"
%include "btllib/data_saveload.hpp"
%include "btllib/status.hpp"
%include "btllib/seq_writer.hpp"
%include "btllib/seq.hpp"

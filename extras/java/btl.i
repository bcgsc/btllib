%module btl

%{
#include "btl/seq_reader.hpp"
#include "btl/graph.hpp"
#include "btl/status.hpp"
#include "btl/seq.hpp"
%}

%include <java.swg>
%include <various.i>
%include <std_string.i>

%include "extra.i"

%include "btl/seq_reader.hpp"
%include "btl/graph.hpp"
%include "btl/status.hpp"
%include "btl/seq.hpp"

%module btl

%{
#include "btl/graph.hpp"
#include "btl/seq_reader.hpp"
#include "btl/seq.hpp"
#include "btl/status.hpp"
%}

%include <java.swg>
%include <various.i>
%include <std_string.i>

%include "extra.i"

%include "btl/graph.hpp"
%include "btl/seq_reader.hpp"
%include "btl/seq.hpp"
%include "btl/status.hpp"

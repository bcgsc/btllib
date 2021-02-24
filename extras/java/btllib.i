%module btllib

%{
#include "btllib/seq_reader.hpp"
#include "btllib/bloom_filter.hpp"
#include "btllib/seq.hpp"
#include "btllib/counting_bloom_filter.hpp"
#include "btllib/indexlr.hpp"
#include "btllib/nthash.hpp"
#include "btllib/seq_writer.hpp"
#include "btllib/data_stream.hpp"
#include "btllib/cstring.hpp"
#include "btllib/order_queue.hpp"
#include "btllib/graph.hpp"
#include "btllib/status.hpp"
#include "btllib/mi_bloom_filter.hpp"
#include "btllib/util.hpp"
%}

%include <java.swg>
%include <various.i>
%include <std_string.i>

%include "extra.i"

%include "btllib/seq_reader.hpp"
%include "btllib/bloom_filter.hpp"
%include "btllib/seq.hpp"
%include "btllib/counting_bloom_filter.hpp"
%include "btllib/indexlr.hpp"
%include "btllib/nthash.hpp"
%include "btllib/seq_writer.hpp"
%include "btllib/data_stream.hpp"
%include "btllib/cstring.hpp"
%include "btllib/order_queue.hpp"
%include "btllib/graph.hpp"
%include "btllib/status.hpp"
%include "btllib/mi_bloom_filter.hpp"
%include "btllib/util.hpp"

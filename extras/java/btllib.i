%module btllib

%{
#include "btllib/graph.hpp"
#include "btllib/seq_writer.hpp"
#include "btllib/data_stream.hpp"
#include "btllib/bloom_filter.hpp"
#include "btllib/seq_reader.hpp"
#include "btllib/util.hpp"
#include "btllib/order_queue.hpp"
#include "btllib/seq.hpp"
#include "btllib/nthash.hpp"
#include "btllib/status.hpp"
#include "btllib/MIBloomFilter.hpp"
%}

%include <java.swg>
%include <various.i>
%include <std_string.i>

%include "extra.i"

%include "btllib/graph.hpp"
%include "btllib/seq_writer.hpp"
%include "btllib/data_stream.hpp"
%include "btllib/bloom_filter.hpp"
%include "btllib/seq_reader.hpp"
%include "btllib/util.hpp"
%include "btllib/order_queue.hpp"
%include "btllib/seq.hpp"
%include "btllib/nthash.hpp"
%include "btllib/status.hpp"
%include "btllib/MIBloomFilter.hpp"

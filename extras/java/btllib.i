%module btllib

%{
#include "btllib/index_queue.hpp"
#include "btllib/seq_reader.hpp"
#include "btllib/graph.hpp"
#include "btllib/util.hpp"
#include "btllib/data_saveload.hpp"
#include "btllib/status.hpp"
#include "btllib/counting_bloom_filter.hpp"
#include "btllib/seq_writer.hpp"
#include "btllib/seq.hpp"
#include "btllib/counting_kmer_set.hpp"
#include "btllib/bloom_filter.hpp"
#include "btllib/kmer_set.hpp"
%}

%include <java.swg>
%include <various.i>
%include <std_string.i>

%include "extra.i"

%include "btllib/index_queue.hpp"
%include "btllib/seq_reader.hpp"
%include "btllib/graph.hpp"
%include "btllib/util.hpp"
%include "btllib/data_saveload.hpp"
%include "btllib/status.hpp"
%include "btllib/counting_bloom_filter.hpp"
%include "btllib/seq_writer.hpp"
%include "btllib/seq.hpp"
%include "btllib/counting_kmer_set.hpp"
%include "btllib/bloom_filter.hpp"
%include "btllib/kmer_set.hpp"

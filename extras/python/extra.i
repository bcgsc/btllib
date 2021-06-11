%rename(__str__) btllib::Graph::to_string;
%feature("python:slot", "tp_str", functype="reprfunc") btllib::Graph::to_string;
%rename(__bool__) btllib::SeqReader::Record::operator bool;
%feature("python:slot", "nb_bool", functype="inquiry") btllib::SeqReader::Record::operator bool;
%rename(__bool__) btllib::Indexlr::Record::operator bool;
%feature("python:slot", "nb_bool", functype="inquiry") btllib::Indexlr::Record::operator bool;

%rename (SeqReaderRecord) btllib::SeqReader::Record;
%rename (SeqReaderFlag) btllib::SeqReader::Flag;
%rename (IndexlrRecord) btllib::Indexlr::Record;
%rename (IndexlrFlag) btllib::Indexlr::Flag;

%rename(__iter__) btllib::SeqReader::begin;
%ignore btllib::SeqReader::end;
%ignore btllib::SeqReader::RecordIterator::operator++;
%ignore btllib::SeqReader::RecordIterator::operator!=;
%ignore btllib::SeqReader::RecordIterator::operator*;
%rename(__next__) btllib::SeqReader::RecordIterator::next;

%feature("python:slot", "tp_iter", functype="getiterfunc") btllib::SeqReader::begin;
%feature("python:slot", "tp_iternext", functype="iternextfunc") btllib::SeqReader::RecordIterator::next;

%extend btllib::SeqReader {
  btllib::SeqReader* __enter__() {
    return $self;
  }
}

%extend btllib::SeqReader {
  void __exit__(PyObject*, PyObject*, PyObject*) {
    $self->close();
  }
}

%exception btllib::SeqReader::RecordIterator::next {
  $action
  if (!bool(result)) {
    PyErr_SetNone(PyExc_StopIteration);
    SWIG_fail;
  }
}

%{
  using btllib::SpacedSeed;
%}

%ignore btllib::CString;
%ignore btllib::CString::operator=;
%ignore btllib::OrderQueue;
%ignore btllib::OrderQueue::Block;
%ignore btllib::OrderQueue::Slot;
%ignore btllib::OrderQueue::Block::operator=;
%ignore btllib::OrderQueue::Block::operator=;
%ignore btllib::OrderQueue::Slot::operator=;
%ignore btllib::OrderQueue::Slot::operator=;

%feature("flatnested", "1");

%ignore btllib::DataStream::operator FILE*() const;
%ignore btllib::DataSource::operator FILE*() const;
%ignore btllib::DataSink::operator FILE*() const;

%ignore btllib::BLOOM_FILTER_MAGIC_HEADER;
%ignore btllib::COUNTING_BLOOM_FILTER_MAGIC_HEADER;

%ignore btllib::SeqReader::read_fasta_buffer;
%ignore btllib::SeqReader::read_fastq_buffer;
%ignore btllib::SeqReader::read_sam_buffer;
%ignore btllib::SeqReader::read_gfa2_buffer;

%ignore btllib::SeqReader::read_fasta_transition;
%ignore btllib::SeqReader::read_fastq_transition;
%ignore btllib::SeqReader::read_sam_transition;
%ignore btllib::SeqReader::read_gfa2_transition;

%ignore btllib::SeqReader::read_fasta_file;
%ignore btllib::SeqReader::read_fastq_file;
%ignore btllib::SeqReader::read_sam_file;
%ignore btllib::SeqReader::read_gfa2_file;

%ignore btllib::NtHash::NtHash(const char*, size_t, unsigned, unsigned, size_t pos = 0);
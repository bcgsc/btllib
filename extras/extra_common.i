%ignore operator<<;

%ignore btllib::DataStream;
%ignore btllib::DataStream::operator FILE*() const;
%ignore btllib::DataSource::operator FILE*() const;
%ignore btllib::DataSink::operator FILE*() const;

%ignore btllib::CString;
%ignore btllib::CString::operator=;

%ignore btllib::OrderQueue;
%ignore btllib::OrderQueue::Block;
%ignore btllib::OrderQueue::Slot;
%ignore btllib::OrderQueue::Block::operator=;
%ignore btllib::OrderQueue::Slot::operator=;

%ignore btllib::SeqReaderFastaModule;
%ignore btllib::SeqReaderFastqModule;
%ignore btllib::SeqReaderGfa2Module;
%ignore btllib::SeqReaderSamModule;
%ignore btllib::SeqReaderMultilineFastaModule;
%ignore btllib::SeqReaderMultilineFastqModule;

%rename (SeqReaderRecord) btllib::SeqReader::Record;
%rename (SeqReaderFlag) btllib::SeqReader::Flag;
%rename (IndexlrRecord) btllib::Indexlr::Record;
%rename (IndexlrFlag) btllib::Indexlr::Flag;

%ignore btllib::IORedirection;
%ignore btllib::IORedirection::in;
%ignore btllib::IORedirection::out;
%ignore btllib::IORedirection::err;

%ignore btllib::ProcessPipeline;
%ignore btllib::ProcessPipeline::in;
%ignore btllib::ProcessPipeline::out;

%ignore btllib::ProcessPipelineInternal;

%ignore btllib::BLOOM_FILTER_MAGIC_HEADER;
%ignore btllib::COUNTING_BLOOM_FILTER_MAGIC_HEADER;

%ignore btllib::SeqReader::RecordIterator::operator++;
%ignore btllib::SeqReader::RecordIterator::operator!=;
%ignore btllib::SeqReader::RecordIterator::operator*;

%ignore btllib::NtHash::NtHash(const char*, size_t, unsigned, unsigned, size_t pos = 0);

%ignore btllib::BloomFilterInitializer;
%ignore btllib::BloomFilterInitializer::operator=;
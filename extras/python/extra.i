%rename(__str__) btllib::Graph::to_string;
%feature("python:slot", "tp_str", functype="reprfunc") btllib::Graph::to_string;
%rename(__bool__) btllib::SeqReader::Record::operator bool;
%feature("python:slot", "nb_bool", functype="inquiry") btllib::SeqReader::Record::operator bool;
%rename(__bool__) btllib::Indexlr::Record::operator bool;
%feature("python:slot", "nb_bool", functype="inquiry") btllib::Indexlr::Record::operator bool;

%rename(__iter__) btllib::SeqReader::begin;
%ignore btllib::SeqReader::RecordIterator::begin;
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

%feature("flatnested", "1");

%ignore btllib::NtHash::NtHash(const std::string&, unsigned, unsigned k, size_t pos = 0);

%template(VectorUint64t) std::vector<uint64_t>;

%{
  static_assert(sizeof(long unsigned int) >= sizeof(uint64_t), "Python wrappers are using wrong size integers.");
%}

%typemap(out) uint64_t* btllib::NtHash::hashes %{
  $result = PyList_New(arg1->get_hash_num());
  for (unsigned i = 0; i < arg1->get_hash_num(); ++i) {
    PyList_SetItem($result, i, PyLong_FromUnsignedLong($1[i]));
  }
%}

%{
  #include <map>
  #include <mutex>
  static long nthash_last_id = 0;
  static std::mutex nthash_mutex;
  static std::map<long, std::string> nthash_strings;
  static std::map<btllib::NtHash*, long> nthash_ids;
%}

%extend btllib::NtHash {
  NtHash(const char* seq, unsigned hash_num, unsigned k, size_t pos = 0)
  {
    std::unique_lock<std::mutex> lock(nthash_mutex);
    nthash_strings[++nthash_last_id] = seq;
    auto *nthash = new btllib::NtHash(nthash_strings[nthash_last_id], hash_num, k);
    nthash_ids[nthash] = nthash_last_id;
    return nthash;
  }

  ~NtHash() {
    std::unique_lock<std::mutex> lock(nthash_mutex);
    nthash_strings.erase(nthash_ids[self]);
    nthash_ids.erase(self);
  }
}
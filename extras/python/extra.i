%rename(__str__) btllib::Graph::to_string;
%feature("python:slot", "tp_str", functype="reprfunc") btllib::Graph::to_string;

%feature("python:slot", "nb_bool", functype="inquiry") btllib::SeqReader::Record::operator bool;
%rename(__bool__) btllib::SeqReader::Record::operator bool;
%feature("python:slot", "nb_bool", functype="inquiry") btllib::Indexlr::Record::operator bool;
%rename(__bool__) btllib::Indexlr::Record::operator bool;

%rename(__iter__) btllib::SeqReader::begin;
%ignore btllib::SeqReader::end;
%rename(__next__) btllib::SeqReader::RecordIterator::next;
%rename(__iter__) btllib::Indexlr::begin;
%ignore btllib::Indexlr::end;
%rename(__next__) btllib::Indexlr::RecordIterator::next;

%feature("python:slot", "tp_iter", functype="getiterfunc") btllib::SeqReader::begin;
%feature("python:slot", "tp_iternext", functype="iternextfunc") btllib::SeqReader::RecordIterator::next;
%feature("python:slot", "tp_iter", functype="getiterfunc") btllib::Indexlr::begin;
%feature("python:slot", "tp_iternext", functype="iternextfunc") btllib::Indexlr::RecordIterator::next;

%extend btllib::SeqReader {
  btllib::SeqReader* __enter__() {
    return $self;
  }

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

%extend btllib::Indexlr {
  btllib::Indexlr* __enter__() {
    return $self;
  }

  void __exit__(PyObject*, PyObject*, PyObject*) {
    $self->close();
  }
}

%exception btllib::Indexlr::RecordIterator::next {
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

%{
  static_assert(sizeof(long unsigned int) >= sizeof(uint64_t), "Python wrappers are using wrong size integers.");
%}

%typemap(out) uint64_t* btllib::NtHash::hashes %{
  $result = PyTuple_New(arg1->get_hash_num());
  for (unsigned i = 0; i < arg1->get_hash_num(); ++i) {
    PyTuple_SetItem($result, i, PyLong_FromUnsignedLong($1[i]));
  }
%}

%typemap(out) uint64_t* btllib::SeedNtHash::hashes %{
  $result = PyTuple_New(arg1->get_hash_num());
  for (unsigned i = 0; i < arg1->get_hash_num(); ++i) {
    PyTuple_SetItem($result, i, PyLong_FromUnsignedLong($1[i]));
  }
%}

%typemap(out) std::vector<btllib::Indexlr::Minimizer>* %{
  $result = PyList_New($1->size());
  for (unsigned i = 0; i < $1->size(); ++i) {
    PyObject *item = SWIG_NewPointerObj(new btllib::Indexlr::Minimizer((*($1))[i]), SWIGTYPE_p_btllib__Indexlr__Minimizer, SWIG_POINTER_OWN);
    PyList_SetItem($result, i, item);
  }
%}
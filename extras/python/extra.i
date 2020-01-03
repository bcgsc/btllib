%feature("python:slot", "tp_str", functype="reprfunc") btl::Graph::to_string;

%ignore btllib::DataSource::operator FILE*();
%ignore btllib::DataSink::operator FILE*();
%ignore btllib::KmerSet::insert(const char*);
%ignore btllib::KmerSet::contains(const char*);

%ignore btllib::CountingKmerSet::insert(const char*);
%ignore btllib::CountingKmerSet::count(const char*);

%ignore btllib::_Data::operator FILE*();
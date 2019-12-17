%ignore operator<<(std::ostream&, const Graph&);
%ignore operator const void*() const;

%feature("python:slot", "tp_str", functype="reprfunc") btl::Graph::to_string;

%inline %{
#include <sstream>
%}

%extend btl::Graph {
    std::string to_string() { std::stringstream ss; ss << *self; return ss.str(); }
}
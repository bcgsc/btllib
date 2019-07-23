%{
#include "sequence.h"
%}

%include <pyprimtypes.swg>
%include <pyopers.swg>
%include <std_common.i>
%include <cstring.i>
%include <std_string.i>
%include <std_iostream.i>

%feature ("flatnested");

%ignore btl::Sequence::Base::operator=(char);
%ignore btl::Sequence::Base::operator=(const Base&);

%ignore btl::Sequence::Base::operator char() const;

%ignore btl::Sequence::Sequence(std::initializer_list<char>);
%ignore btl::Sequence::Sequence(Sequence&&);

%ignore btl::Sequence::operator=(const Sequence&);
%ignore btl::Sequence::operator=(std::string);
%ignore btl::Sequence::operator=(const char*);
%ignore btl::Sequence::operator=(std::initializer_list<char>);
%ignore btl::Sequence::operator=(Sequence&&);
%ignore btl::Sequence::operator=(std::string&&);

%ignore btl::Sequence::operator const std::string&() const noexcept;
%ignore btl::Sequence::operator const char*() const noexcept;

%ignore btl::Sequence::operator+=(std::initializer_list<char>);

%ignore btl::Sequence::operator[](size_t pos);
%ignore btl::Sequence::operator[](size_t pos) const;
%ignore btl::Sequence::at(size_t pos) const;
%ignore btl::Sequence::back() const;
%ignore btl::Sequence::front() const;

%ignore btl::Sequence::append(std::initializer_list<char>);
%ignore btl::Sequence::assign(std::initializer_list<char>);
%ignore btl::Sequence::assign(Sequence&&) noexcept;
%ignore btl::Sequence::assign(std::string&&) noexcept;
%ignore btl::Sequence::insert(const_iterator, std::initializer_list<char>);
%ignore btl::Sequence::replace(const_iterator, const_iterator, std::initializer_list<char>);

%ignore btl::swap(Sequence&, Sequence&);

%ignore operator+(const Sequence&, const Sequence&);
%ignore operator+(Sequence&&, Sequence&&);
%ignore operator+(Sequence&&, const Sequence&);
%ignore operator+(const Sequence&, Sequence&&);
%ignore operator+(const std::string&, const Sequence&);
%ignore operator+(const Sequence&, const std::string&);
%ignore operator+(const std::string&, Sequence&&);
%ignore operator+(const Sequence&, std::string&&);
%ignore operator+(std::string&&, const Sequence&);
%ignore operator+(Sequence&&, const std::string&);
%ignore operator+(Sequence&&, std::string&&);
%ignore operator+(std::string&&, Sequence&&);
%ignore operator+(const char*, const Sequence&);
%ignore operator+(const Sequence&, const char*);
%ignore operator+(const char*, Sequence&&);
%ignore operator+(Sequence&&, const char*);
%ignore operator+(char, const Sequence&);
%ignore operator+(const Sequence&, char);
%ignore operator+(char, Sequence&&);
%ignore operator+(Sequence&&, char);

%ignore btl::operator+(const Sequence&, const Sequence&);
%ignore btl::operator+(Sequence&&, Sequence&&);
%ignore btl::operator+(Sequence&&, Sequence&&);
%ignore btl::operator+(Sequence&&, const Sequence&);
%ignore btl::operator+(const Sequence&, Sequence&&);
%ignore btl::operator+(const std::string&, const Sequence&);
%ignore btl::operator+(const Sequence&, const std::string&);
%ignore btl::operator+(const std::string&, Sequence&&);
%ignore btl::operator+(const Sequence&, std::string&&);
%ignore btl::operator+(std::string&&, const Sequence&);
%ignore btl::operator+(Sequence&&, const std::string&);
%ignore btl::operator+(Sequence&&, std::string&&);
%ignore btl::operator+(std::string&&, Sequence&&);
%ignore btl::operator+(const char*, const Sequence&);
%ignore btl::operator+(const Sequence&, const char*);
%ignore btl::operator+(const char*, Sequence&&);
%ignore btl::operator+(Sequence&&, const char*);
%ignore btl::operator+(char, const Sequence&);
%ignore btl::operator+(const Sequence&, char);
%ignore btl::operator+(char, Sequence&&);
%ignore btl::operator+(Sequence&&, char);

%ignore operator==(const Sequence&, const Sequence&);
%ignore operator==(const std::string&, const Sequence&);
%ignore operator==(const Sequence&, const std::string&);
%ignore operator==(const char*, const Sequence&);
%ignore operator==(const Sequence&, const char*);

%ignore operator!=(const Sequence&, const Sequence&);
%ignore operator!=(const std::string&, const Sequence&);
%ignore operator!=(const Sequence&, const std::string&);
%ignore operator!=(const char*, const Sequence&);
%ignore operator!=(const Sequence&, const char*);

%ignore operator<(const Sequence&, const Sequence&);
%ignore operator<(const std::string&, const Sequence&);
%ignore operator<(const Sequence&, const std::string&);
%ignore operator<(const char*, const Sequence&);
%ignore operator<(const Sequence&, const char*);

%ignore operator<=(const Sequence&, const Sequence&);
%ignore operator<=(const std::string&, const Sequence&);
%ignore operator<=(const Sequence&, const std::string&);
%ignore operator<=(const char*, const Sequence&);
%ignore operator<=(const Sequence&, const char*);

%ignore operator>(const Sequence&, const Sequence&);
%ignore operator>(const std::string&, const Sequence&);
%ignore operator>(const Sequence&, const std::string&);
%ignore operator>(const char*, const Sequence&);
%ignore operator>(const Sequence&, const char*);

%ignore operator>=(const Sequence&, const Sequence&);
%ignore operator>=(const std::string&, const Sequence&);
%ignore operator>=(const Sequence&, const std::string&);
%ignore operator>=(const char*, const Sequence&);
%ignore operator>=(const Sequence&, const char*);

%ignore operator>>(std::istream&, Sequence&);
%ignore operator<<(std::ostream&, const Sequence&);

%ignore getline(std::istream&, Sequence&, char);
%ignore getline(std::istream&&, Sequence&, char);
%ignore getline(std::istream&, Sequence&);
%ignore getline(std::istream&&, Sequence&);

%extend btl::Sequence {
    char __getitem__(size_t key) { return self->operator[](key); }
    void __setitem__(size_t key, char value) { self->operator[](key) = value; }

    std::string __str__() { return std::string(*self); }

    Sequence __add__(const Sequence& seq) { return *self + seq; }
    Sequence __radd__(const Sequence& seq) { return seq + *self; }
    Sequence __add__(const std::string& seq) { return *self + seq; }
    Sequence __radd__(const std::string& seq) { return seq + *self; }
    Sequence __add__(const char* seq) { return *self + seq; }
    Sequence __radd__(const char* seq) { return seq + *self; }
    Sequence __add__(char base) { return *self + base; }
    Sequence __radd__(char base) { return base + *self; }
    bool __eq__(const Sequence& seq) { return *self == seq; }
    bool __eq__(const std::string& seq) { return *self == seq; }
    bool __eq__(const char* seq) { return *self == seq; }
    bool __ne__(const Sequence& seq) { return *self != seq; }
    bool __ne__(const std::string& seq) { return *self != seq; }
    bool __ne__(const char* seq) { return *self != seq; }
    bool __gt__(const Sequence& seq) { return *self > seq; }
    bool __gt__(const std::string& seq) { return *self > seq; }
    bool __gt__(const char* seq) { return *self > seq; }
    bool __ge__(const Sequence& seq) { return *self >= seq; }
    bool __ge__(const std::string& seq) { return *self >= seq; }
    bool __ge__(const char* seq) { return *self >= seq; }
    bool __lt__(const Sequence& seq) { return *self < seq; }
    bool __lt__(const std::string& seq) { return *self < seq; }
    bool __lt__(const char* seq) { return *self < seq; }
    bool __le__(const Sequence& seq) { return *self <= seq; }
    bool __le__(const std::string& seq) { return *self <= seq; }
    bool __le__(const char* seq) { return *self <= seq; }
};

%include "sequence.h"

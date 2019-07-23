%{
#define SWIG_FILE_WITH_INIT
#include "sequence.h"
%}

%include <java.swg>
%include <various.i>
%include <std_string.i>

%feature ("flatnested");

%ignore btl::Sequence::Base::operator=(char);
%ignore btl::Sequence::Base::operator=(const Base&);

%ignore btl::Sequence::Base::operator char() const;

%ignore btl::Sequence::Base::operator~() const;

%ignore btl::Sequence::Sequence(const char*);
%ignore btl::Sequence::Sequence(const char*, size_t);
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
%rename(add) btl::Sequence::operator+=(const Sequence&);
%rename(add) btl::Sequence::operator+=(const std::string&);
%ignore btl::Sequence::operator+=(const char*);
%rename(add) btl::Sequence::operator+=(char);

%ignore btl::Sequence::begin() const noexcept;
%ignore btl::Sequence::end() const noexcept;
%ignore btl::Sequence::rbegin() const noexcept;
%ignore btl::Sequence::rend() const noexcept;

%ignore btl::Sequence::operator[](size_t pos);
%ignore btl::Sequence::operator[](size_t pos) const;
%ignore btl::Sequence::at(size_t pos) const;
%ignore btl::Sequence::back() const;
%ignore btl::Sequence::front() const;

%ignore btl::Sequence::append(const char*);
%ignore btl::Sequence::append(std::initializer_list<char>);
%ignore btl::Sequence::assign(const char*);
%ignore btl::Sequence::assign(const char*, size_t);
%ignore btl::Sequence::assign(std::initializer_list<char>);
%ignore btl::Sequence::assign(Sequence&&) noexcept;
%ignore btl::Sequence::assign(std::string&&) noexcept;
%ignore btl::Sequence::insert(size_t, const char*);
%ignore btl::Sequence::insert(size_t, const char*, size_t);
%ignore btl::Sequence::insert(const_iterator, std::initializer_list<char>);
%ignore btl::Sequence::replace(size_t, size_t, const char*);
%ignore btl::Sequence::replace(const_iterator, const_iterator, const char*);
%ignore btl::Sequence::replace(size_t, size_t, const char*, size_t);
%ignore btl::Sequence::replace(const_iterator, const_iterator, std::initializer_list<char>);

%ignore btl::swap(Sequence&, Sequence&);

%ignore btl::Sequence::find(const char*, size_t) const;
%ignore btl::Sequence::find(const char*) const;
%ignore btl::Sequence::rfind(const char*, size_t) const;
%ignore btl::Sequence::rfind(const char*) const;
%ignore btl::Sequence::find_first_of(const char*, size_t) const;
%ignore btl::Sequence::find_first_of(const char*) const;
%ignore btl::Sequence::find_last_of(const char*, size_t) const;
%ignore btl::Sequence::find_last_of(const char*) const;
%ignore btl::Sequence::find_first_not_of(const char*, size_t) const;
%ignore btl::Sequence::find_first_not_of(const char*) const;
%ignore btl::Sequence::find_last_not_of(const char*, size_t) const;
%ignore btl::Sequence::find_last_not_of(const char*) const;

%ignore btl::Sequence::compare(const char*) const;
%ignore btl::Sequence::compare(size_t, size_t, const char*) const;
%ignore btl::Sequence::compare(size_t, size_t, const char*, size_t) const;

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

%include "sequence.h"

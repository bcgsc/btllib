#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "base.h"

#include <string>
#include <algorithm>
#include <cassert>
#include <iostream>

namespace btl {

class Sequence {

public:

    Sequence();
    Sequence(std::string bases);
    Sequence(const char* bases);
    Sequence(const Sequence& seq);
    Sequence(Sequence&& seq);

    Sequence& operator=(std::string bases);
    Sequence& operator=(const char* bases);
    Sequence& operator=(const Sequence& seq);
    Sequence& operator=(Sequence&& seq);

    void reverseComplement();

    Sequence getReverseComplement();

    Sequence operator~();

    static inline const size_t& npos = std::string::npos;

    using iterator = std::string::iterator;
    using const_iterator = std::string::const_iterator;
    using reverse_iterator = std::string::reverse_iterator;
    using const_reverse_iterator = std::string::const_reverse_iterator;

    iterator begin() noexcept;
    const_iterator begin() const noexcept;
    iterator end() noexcept;
    const_iterator end() const noexcept;
    reverse_iterator rbegin() noexcept;
    const_reverse_iterator rbegin() const noexcept;
    reverse_iterator rend() noexcept;
    const_reverse_iterator rend() const noexcept;
    const_iterator cbegin() const noexcept;
    const_iterator cend() const noexcept;
    const_reverse_iterator crbegin() const noexcept;
    const_reverse_iterator crend() const noexcept;

    size_t size() const noexcept;
    size_t length() const noexcept;
    size_t max_size() const noexcept;
    void resize(size_t n);
    void resize(size_t n, Base base);
    size_t capacity() const noexcept;
    void reserve(size_t n = 0);
    void clear() noexcept;
    bool empty() const noexcept;
    void shrink_to_fit();

    Base& operator[](size_t pos);
    const Base& operator[](size_t pos) const;
    Base& at(size_t pos);
    const Base& at(size_t pos) const;
    Base& back();
    const Base& back() const;
    Base& front();
    const Base& front() const;

    Sequence& operator+=(const Sequence& rhs);
    Sequence& operator+=(Base rhs);
    Sequence& operator+=(std::initializer_list<char> rhs);

    Sequence& append(const Sequence& seq);
    Sequence& append(const Sequence& seq, size_t subpos, size_t sublen);
    Sequence& append(const std::string& bases);
    Sequence& append(const std::string& bases, size_t subpos, size_t sublen);
    Sequence& append(const char* bases);
    Sequence& append(const char* bases, size_t n);
    Sequence& append(size_t n, Base base);
    template <class InputIterator>
    Sequence& append(InputIterator first, InputIterator last);
    Sequence& append(std::initializer_list<char> bases);
    void push_back (char c);
    Sequence& assign(const Sequence& seq);
    Sequence& assign(const Sequence& seq, size_t subpos, size_t sublen = npos);
    Sequence& assign(const std::string& bases);
    Sequence& assign(const std::string& bases, size_t subpos, size_t sublen = npos);
    Sequence& assign(const char* bases);
    Sequence& assign(const char* bases, size_t n);
    Sequence& assign(size_t n, char c);
    template <class InputIterator>
    Sequence& assign(InputIterator first, InputIterator last);
    Sequence& assign(std::initializer_list<char> bases);
    Sequence& assign(Sequence&& seq) noexcept;
    Sequence& assign(std::string&& bases) noexcept;

    string& insert (size_t pos, const string& str);
    string& insert (size_t pos, const string& str, size_t subpos, size_t sublen = npos);
    string& insert (size_t pos, const char* s);
    string& insert (size_t pos, const char* s, size_t n);
    string& insert (size_t pos,   size_t n, char c);
    iterator insert (const_iterator p, size_t n, char c);
    iterator insert (const_iterator p, char c);
    template <class InputIterator>
    iterator insert (iterator p, InputIterator first, InputIterator last);
    string& insert (const_iterator p, initializer_list<char> il);

    friend Sequence operator+(const Sequence& lhs, const Sequence& rhs);
    friend Sequence operator+(const std::string& lhs, const Sequence& rhs);
    friend Sequence operator+(const Sequence& lhs, const std::string& rhs);
    friend Sequence operator+(const char* lhs, const Sequence& rhs);
    friend Sequence operator+(const Sequence& lhs, const char* rhs);
    friend Sequence operator+(char lhs, const Sequence& rhs);
    friend Sequence operator+(const Sequence& lhs, char rhs);

    friend bool operator==(const Sequence& lhs, const Sequence& rhs);
    friend bool operator==(const std::string& lhs, const Sequence& rhs);
    friend bool operator==(const Sequence& lhs, const std::string& rhs);
    friend bool operator==(const char* lhs, const Sequence& rhs);
    friend bool operator==(const Sequence& lhs, const char* rhs);
    friend bool operator==(char lhs, const Sequence& rhs);
    friend bool operator==(const Sequence& lhs, char rhs);

    friend bool operator!=(const Sequence& lhs, const Sequence& rhs);
    friend bool operator!=(const std::string& lhs, const Sequence& rhs);
    friend bool operator!=(const Sequence& lhs, const std::string& rhs);
    friend bool operator!=(const char* lhs, const Sequence& rhs);
    friend bool operator!=(const Sequence& lhs, const char* rhs);
    friend bool operator!=(char lhs, const Sequence& rhs);
    friend bool operator!=(const Sequence& lhs, char rhs);

    friend std::istream& operator>>(std::istream& is, Sequence& rhs);
    friend std::ostream& operator<<(std::ostream& os, const Sequence& rhs);

    static void check_validity(const std::string& bases);
    static void check_validity(const char* bases);
    static void check_validity(std::initializer_list<char> bases);

private:

    void set(std::string bases);
    void set(const char* bases);

    std::string s;

};

inline Sequence::Sequence() {}
inline Sequence::Sequence(std::string bases): s(std::move(bases)) { check_validity(s); }
inline Sequence::Sequence(const char* bases): s(bases) { check_validity(s); }
inline Sequence::Sequence(const Sequence& seq): s(seq.s) {}
inline Sequence::Sequence(Sequence&& seq): s(seq.s) {}

inline Sequence& Sequence::operator=(std::string bases) { s = std::move(bases); check_validity(s); return *this; }
inline Sequence& Sequence::operator=(const char* bases) { s = std::string(bases); check_validity(s); return *this; }
inline Sequence& Sequence::operator=(const Sequence& seq) { s = seq.s; return *this; }
inline Sequence& Sequence::operator=(Sequence&& seq) { s = seq.s; return *this; }

inline void Sequence::reverseComplement() {
	std::reverse(s.begin(), s.end());
	std::transform(s.begin(), s.end(), s.begin(),
        [] (unsigned char c) { return COMPLEMENTS[c]; }
    );
}

inline Sequence Sequence::getReverseComplement() {
    Sequence seq;
    seq.reserve(s.size());
    std::for_each(s.rbegin(), s.rend(),
        [&] (char c) { seq.s += COMPLEMENTS[(unsigned char)c]; }
    );
    return seq;
}

inline Sequence Sequence::operator~() {
    return getReverseComplement();
}

inline Sequence::iterator Sequence::begin() noexcept { return s.begin(); }
inline Sequence::const_iterator Sequence::begin() const noexcept { return s.begin(); }
inline Sequence::iterator Sequence::end() noexcept { return s.end(); }
inline Sequence::const_iterator Sequence::end() const noexcept { return s.end(); }
inline Sequence::reverse_iterator Sequence::rbegin() noexcept { return s.rbegin(); }
inline Sequence::const_reverse_iterator Sequence::rbegin() const noexcept { return s.rbegin(); }
inline Sequence::reverse_iterator Sequence::rend() noexcept { return s.rend(); }
inline Sequence::const_reverse_iterator Sequence::rend() const noexcept { return s.rend(); }
inline Sequence::const_iterator Sequence::cbegin() const noexcept { return s.cbegin(); }
inline Sequence::const_iterator Sequence::cend() const noexcept { return s.cend(); }
inline Sequence::const_reverse_iterator Sequence::crbegin() const noexcept { return s.crbegin(); }
inline Sequence::const_reverse_iterator Sequence::crend() const noexcept { return s.crend(); }

inline size_t Sequence::size() const noexcept { return s.size(); }
inline size_t Sequence::length() const noexcept { return s.length(); }
inline size_t Sequence::max_size() const noexcept { return s.max_size(); }
inline void Sequence::resize (size_t n) { s.resize(n); }
inline void Sequence::resize (size_t n, Base b) { s.resize(n, b); }
inline size_t Sequence::capacity() const noexcept { return s.capacity(); }
inline void Sequence::reserve(size_t n) { s.reserve(n); }
inline void Sequence::clear() noexcept { s.clear(); }
inline bool Sequence::empty() const noexcept { return s.empty(); }
inline void Sequence::shrink_to_fit() { s.shrink_to_fit(); }

inline char& Sequence::operator[](size_t pos) { return s[pos]; }
inline const char& Sequence::operator[](size_t pos) const { return s[pos]; }
inline char& Sequence::at(size_t pos) { return s.at(pos); }
inline const char& Sequence::at(size_t pos) const { return s.at(pos); }
inline char& Sequence::back() { return s.back(); }
inline const char& Sequence::back() const { return s.back(); }
inline char& Sequence::front() { return s.front(); }
inline const char& Sequence::front() const { return s.front(); }

inline Sequence& Sequence::operator+=(const Sequence& rhs) { s += rhs.s; return *this; }
inline Sequence& Sequence::operator+=(const std::string& rhs) { check_validity(rhs); s += rhs; return *this; }
inline Sequence& Sequence::operator+=(const char* rhs) { check_validity(rhs); s += rhs; return *this; }
inline Sequence& Sequence::operator+=(char rhs) { check_validity(rhs); s += rhs; return *this; }
inline Sequence& Sequence::operator+=(std::initializer_list<char> rhs) { check_validity(rhs); s += rhs; return *this; }

inline Sequence& Sequence::append(const Sequence& seq) { s.append(seq.s); return *this; }
inline Sequence& Sequence::append(const Sequence& seq, size_t subpos, size_t sublen) { s.append(seq.s, subpos, sublen); return *this; }
inline Sequence& Sequence::append(const std::string& bases) { check_validity(bases); s.append(bases); return *this; }
inline Sequence& Sequence::append(const std::string& bases, size_t subpos, size_t sublen) { check_validity(bases); s.append(bases, subpos, sublen); return *this; }
inline Sequence& Sequence::append(const char* bases) { check_validity(bases); s.append(bases); return *this; }
inline Sequence& Sequence::append(const char* bases, size_t n) { check_validity(bases); s.append(bases, n); return *this; }
inline Sequence& Sequence::append(size_t n, char c) { check_validity(c); s.append(n, c); return *this; }
template <class InputIterator>
inline Sequence& Sequence::append(InputIterator first, InputIterator last) { s.append(first, last); check_validity(s); return *this; }
inline Sequence& Sequence::append(std::initializer_list<char> bases) { check_validity(bases); s.append(bases); return *this; }
inline void Sequence::push_back(char c) { check_validity(c); s.push_back(c); }
inline Sequence& Sequence::assign(const Sequence& seq) { s.assign(seq.s); return *this; }
inline Sequence& Sequence::assign(const Sequence& seq, size_t subpos, size_t sublen = npos) { s.assign(seq.s, subpos, sublen); return *this; }
inline Sequence& Sequence::assign(const std::string& bases) { check_validity(bases); s.assign(bases); return *this; }
inline Sequence& Sequence::assign(const std::string& bases, size_t subpos, size_t sublen = npos) { check_validity(bases); s.assign(bases, subpos, npos); return *this; }
inline Sequence& Sequence::assign(const char* bases) { check_validity(bases); s.assign(bases); return *this; }
inline Sequence& Sequence::assign(const char* bases, size_t n) { check_validity(bases); s.assign(bases, n); return *this; }
inline Sequence& Sequence::assign(size_t n, char c) { check_validity(c); s.assign(n, c); return *this; }
template <class InputIterator>
inline Sequence& Sequence::assign(InputIterator first, InputIterator last) { s.assign(first, last); check_validity(s); return *this; }
inline Sequence& Sequence::assign(std::initializer_list<char> bases) { check_validity(bases); s.assign(bases); return *this; }
inline Sequence& Sequence::assign(Sequence&& seq) noexcept { s.assign(seq.s); }
inline Sequence& Sequence::assign(std::string&& bases) noexcept { check_validity(bases); s.assign(bases); }

inline Sequence operator+(const Sequence& lhs, const Sequence& rhs) { Sequence seq; seq.set(lhs.s + rhs.s); return seq; }
inline Sequence operator+(const std::string& lhs, const Sequence& rhs) { Sequence::check_validity(lhs); Sequence seq; seq.set(lhs + rhs.s); return seq; }
inline Sequence operator+(const Sequence& lhs, const std::string& rhs) { Sequence::check_validity(rhs); Sequence seq; seq.set(lhs.s + rhs); return seq; }
inline Sequence operator+(const char* lhs, const Sequence& rhs) { Sequence::check_validity(lhs); Sequence seq; seq.set(lhs + rhs.s); return seq; }
inline Sequence operator+(const Sequence& lhs, const char* rhs) { Sequence::check_validity(rhs); Sequence seq; seq.set(lhs.s + rhs); return seq; }
inline Sequence operator+(char lhs, const Sequence& rhs) { Sequence::check_validity(lhs); Sequence seq; seq.set(lhs + rhs.s); return seq; }
inline Sequence operator+(const Sequence& lhs, char rhs) {  Sequence::check_validity(rhs); Sequence seq; seq.set(lhs.s + rhs); return seq; }

inline bool operator==(const Sequence& lhs, const Sequence& rhs) { return lhs.s == rhs.s; }
inline bool operator==(const std::string& lhs, const Sequence& rhs) { return lhs == rhs.s; }
inline bool operator==(const Sequence& lhs, const std::string& rhs) { return lhs.s == rhs; }
inline bool operator==(const char* lhs, const Sequence& rhs) { return lhs == rhs.s; }
inline bool operator==(const Sequence& lhs, const char* rhs) { return lhs.s == rhs; }
inline bool operator==(char lhs, const Sequence& rhs) { return lhs == rhs.s; }
inline bool operator==(const Sequence& lhs, char rhs) { return lhs.s == rhs; }

inline bool operator!=(const Sequence& lhs, const Sequence& rhs) { return lhs.s != rhs.s; }
inline bool operator!=(const std::string& lhs, const Sequence& rhs) { return lhs != rhs.s; }
inline bool operator!=(const Sequence& lhs, const std::string& rhs) { return lhs.s != rhs; }
inline bool operator!=(const char* lhs, const Sequence& rhs) { return lhs != rhs.s; }
inline bool operator!=(const Sequence& lhs, const char* rhs) { return lhs.s != rhs; }
inline bool operator!=(char lhs, const Sequence& rhs) { return lhs != rhs.s; }
inline bool operator!=(const Sequence& lhs, char rhs) { return lhs.s != rhs; }

inline std::istream& operator>>(std::istream& is, Sequence& rhs) {
    is >> rhs.s;
    return is;
}

inline std::ostream& operator<<(std::ostream& os, const Sequence& rhs) {
    os << rhs.s;
    return os;
}

inline void Sequence::set(std::string bases) { s = std::move(bases); }
inline void Sequence::set(const char* bases) { s = std::string(bases); }

inline void Sequence::check_validity(const std::string& bases) {
    std::for_each(bases.begin(), bases.end(),
        [] (char c) { Base::check_validity(c); }
    );
}

inline void Sequence::check_validity(const char* bases) {
    for (unsigned i = 0; bases[i] != '\0'; ++i) {
        assert(Base::check_validity(bases[i]););   
    }
}

inline void Sequence::check_validity(std::initializer_list<char> bases) {
    std::for_each(bases.begin(), bases.end(),
        [] (char c) { Base::check_validity(c); }
    );
}

}

#endif
#ifndef BASIC_SEQUENCE_H
#define BASIC_SEQUENCE_H

#include <string>
#include <algorithm>
#include <cassert>
#include <iostream>

namespace btl {

class BasicSequence {

public:

    class Base {

    public:

        Base& operator=(char base);
        Base& operator=(Base base);

        operator char() const;

        void complement();
        char operator~() const;

        static inline void validate(char base);

        friend class BasicSequence;

        Base(char& base);
        Base(Base& base);

        char& b;
        static const char COMPLEMENTS[256];
    };

    static inline const size_t& npos = std::string::npos;

    BasicSequence();
    BasicSequence(const BasicSequence& seq);
    BasicSequence(std::string bases);
    BasicSequence(const BasicSequence& seq, size_t pos, size_t len = npos);
    BasicSequence(const std::string& bases, size_t pos, size_t len = npos);
    BasicSequence(const char* bases);
    BasicSequence(const char* bases, size_t n);
    BasicSequence(size_t n, char base);
    template <class InputIterator>
    BasicSequence(InputIterator first, InputIterator last);
    BasicSequence(std::initializer_list<char> bases);
    BasicSequence(BasicSequence&& seq);

    BasicSequence& operator=(const BasicSequence& seq);
    BasicSequence& operator=(std::string bases);
    BasicSequence& operator=(const char* bases);
    BasicSequence& operator=(std::initializer_list<char> bases);
    BasicSequence& operator=(BasicSequence&& seq);
    BasicSequence& operator=(std::string&& bases);

    operator const std::string&() const noexcept;
    operator const char*() const noexcept;

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
    void resize(size_t n, char base);
    size_t capacity() const noexcept;
    void reserve(size_t n = 0);
    void clear() noexcept;
    bool empty() const noexcept;
    void shrink_to_fit();

    Base operator[](size_t pos);
    const char& operator[](size_t pos) const;
    Base at(size_t pos);
    const char& at(size_t pos) const;
    Base back();
    const char& back() const;
    Base front();
    const char& front() const;

    BasicSequence& operator+=(const BasicSequence& rhs);
    BasicSequence& operator+=(const std::string& rhs);
    BasicSequence& operator+=(const char* rhs);
    BasicSequence& operator+=(char rhs);
    BasicSequence& operator+=(std::initializer_list<char> rhs);

    BasicSequence& append(const BasicSequence& seq);
    BasicSequence& append(const BasicSequence& seq, size_t subpos, size_t sublen);
    BasicSequence& append(const std::string& bases);
    BasicSequence& append(const std::string& bases, size_t subpos, size_t sublen);
    BasicSequence& append(const char* bases);
    BasicSequence& append(const char* bases, size_t n);
    BasicSequence& append(size_t n, char base);
    template <class InputIterator>
    BasicSequence& append(InputIterator first, InputIterator last);
    BasicSequence& append(std::initializer_list<char> bases);

    void push_back(char base);

    BasicSequence& assign(const BasicSequence& seq);
    BasicSequence& assign(const BasicSequence& seq, size_t subpos, size_t sublen = npos);
    BasicSequence& assign(const std::string& bases);
    BasicSequence& assign(const std::string& bases, size_t subpos, size_t sublen = npos);
    BasicSequence& assign(const char* bases);
    BasicSequence& assign(const char* bases, size_t n);
    BasicSequence& assign(size_t n, char base);
    template <class InputIterator>
    BasicSequence& assign(InputIterator first, InputIterator last);
    BasicSequence& assign(std::initializer_list<char> bases);
    BasicSequence& assign(BasicSequence&& seq) noexcept;
    BasicSequence& assign(std::string&& bases) noexcept;

    BasicSequence& insert(size_t pos, const BasicSequence& seq);
    BasicSequence& insert(size_t pos, const BasicSequence& seq, size_t subpos, size_t sublen = npos);
    BasicSequence& insert(size_t pos, const std::string& bases);
    BasicSequence& insert(size_t pos, const std::string& bases, size_t subpos, size_t sublen = npos);
    BasicSequence& insert(size_t pos, const char* bases);
    BasicSequence& insert(size_t pos, const char* bases, size_t n);
    BasicSequence& insert(size_t pos, size_t n, char base);
    iterator insert(const_iterator p, size_t n, char base);
    iterator insert(const_iterator p, char base);
    template <class InputIterator>
    iterator insert(iterator p, InputIterator first, InputIterator last);
    BasicSequence& insert(const_iterator p, std::initializer_list<char> bases);

    BasicSequence& erase(size_t pos = 0, size_t len = npos);
    iterator erase(iterator p);
    iterator erase(iterator first, iterator last);

    BasicSequence& replace(size_t pos, size_t len, const BasicSequence& seq);
    BasicSequence& replace(size_t pos, size_t len, const std::string& bases);
    BasicSequence& replace(const_iterator i1, const_iterator i2, const BasicSequence& seq);
    BasicSequence& replace(const_iterator i1, const_iterator i2, const std::string& bases);
    BasicSequence& replace(size_t pos, size_t len, const BasicSequence& seq, size_t subpos, size_t sublen = npos);
    BasicSequence& replace(size_t pos, size_t len, const std::string& bases, size_t subpos, size_t sublen = npos);
    BasicSequence& replace(size_t pos, size_t len, const char* bases);
    BasicSequence& replace(const_iterator i1, const_iterator i2, const char* bases);
    BasicSequence& replace(size_t pos, size_t len, const char* bases, size_t n);
    BasicSequence& replace(const_iterator i1, const_iterator i2, const char* bases, size_t n);
    BasicSequence& replace(size_t pos, size_t len, size_t n, char base);
    BasicSequence& replace(const_iterator i1, const_iterator i2, size_t n, char base);
    template <class InputIterator>
    BasicSequence& replace(const_iterator i1, const_iterator i2, InputIterator first, InputIterator last);
    BasicSequence& replace(const_iterator i1, const_iterator i2, std::initializer_list<char> bases);

    void swap(BasicSequence& seq);
    void swap(std::string& bases);

    void pop_back();

    const char* c_str() const noexcept;
    const char* data() const noexcept;

    using allocator_type = std::string::allocator_type;

    allocator_type get_allocator() const noexcept;

    size_t copy(char* bases, size_t len, size_t pos = 0) const;

    using size_type = std::string::size_type;

    size_t find(const BasicSequence& seq, size_t pos = 0) const noexcept;
    size_t find(const std::string& bases, size_t pos = 0) const noexcept;
    size_t find(const char* bases, size_t pos = 0) const;
    size_t find(const char* bases, size_t pos, size_type n) const;
    size_t find(char base, size_t pos = 0) const noexcept;

    size_t rfind(const BasicSequence& seq, size_t pos = npos) const noexcept;
    size_t rfind(const std::string& bases, size_t pos = npos) const noexcept;
    size_t rfind(const char* bases, size_t pos = npos) const;
    size_t rfind(const char* bases, size_t pos, size_t n) const;
    size_t rfind(char base, size_t pos = npos) const noexcept;

    size_t find_first_of(const BasicSequence& seq, size_t pos = 0) const noexcept;
    size_t find_first_of(const std::string& bases, size_t pos = 0) const noexcept;
    size_t find_first_of(const char* bases, size_t pos = 0) const;
    size_t find_first_of(const char* bases, size_t pos, size_t n) const;
    size_t find_first_of(char base, size_t pos = 0) const noexcept;

    size_t find_last_of(const BasicSequence& seq, size_t pos = npos) const noexcept;
    size_t find_last_of(const std::string& bases, size_t pos = npos) const noexcept;
    size_t find_last_of(const char* bases, size_t pos = npos) const;
    size_t find_last_of(const char* bases, size_t pos, size_t n) const;
    size_t find_last_of(char base, size_t pos = npos) const noexcept;

    size_t find_first_not_of(const BasicSequence& seq, size_t pos = 0) const noexcept;
    size_t find_first_not_of(const std::string& bases, size_t pos = 0) const noexcept;
    size_t find_first_not_of(const char* bases, size_t pos = 0) const;
    size_t find_first_not_of(const char* bases, size_t pos, size_t n) const;
    size_t find_first_not_of(char base, size_t pos = 0) const noexcept;

    size_t find_last_not_of(const BasicSequence& seq, size_t pos = npos) const noexcept;
    size_t find_last_not_of(const std::string& bases, size_t pos = npos) const noexcept;
    size_t find_last_not_of(const char* bases, size_t pos = npos) const;
    size_t find_last_not_of(const char* bases, size_t pos, size_t n) const;
    size_t find_last_not_of(char base, size_t pos = npos) const noexcept;

    BasicSequence substr(size_t pos = 0, size_t len = npos) const;

    int compare(const BasicSequence& seq) const noexcept;
    int compare(const std::string& bases) const noexcept;
    int compare(size_t pos, size_t len, const BasicSequence& seq) const;
    int compare(size_t pos, size_t len, const std::string& bases) const;
    int compare(size_t pos, size_t len, const BasicSequence& seq, size_t subpos, size_t sublen = npos) const;
    int compare(size_t pos, size_t len, const std::string& bases, size_t subpos, size_t sublen = npos) const;
    int compare(const char* bases) const;
    int compare(size_t pos, size_t len, const char* bases) const;
    int compare(size_t pos, size_t len, const char* bases, size_t n) const;

    friend BasicSequence operator+(const BasicSequence& lhs, const BasicSequence& rhs);
    friend BasicSequence operator+(BasicSequence&& lhs, BasicSequence&& rhs);
    friend BasicSequence operator+(BasicSequence&& lhs, const BasicSequence& rhs);
    friend BasicSequence operator+(const BasicSequence& lhs, BasicSequence&& rhs);
    friend BasicSequence operator+(const std::string& lhs, const BasicSequence& rhs);
    friend BasicSequence operator+(const BasicSequence& lhs, const std::string& rhs);
    friend BasicSequence operator+(const std::string& lhs, BasicSequence&& rhs);
    friend BasicSequence operator+(const BasicSequence& lhs, std::string&& rhs);
    friend BasicSequence operator+(std::string&& lhs, const BasicSequence& rhs);
    friend BasicSequence operator+(BasicSequence&& lhs, const std::string& rhs);
    friend BasicSequence operator+(BasicSequence&& lhs, std::string&& rhs);
    friend BasicSequence operator+(std::string&& lhs, BasicSequence&& rhs);
    friend BasicSequence operator+(const char* lhs, const BasicSequence& rhs);
    friend BasicSequence operator+(const BasicSequence& lhs, const char* rhs);
    friend BasicSequence operator+(const char* lhs, BasicSequence&& rhs);
    friend BasicSequence operator+(BasicSequence&& lhs, const char* rhs);
    friend BasicSequence operator+(char lhs, const BasicSequence& rhs);
    friend BasicSequence operator+(const BasicSequence& lhs, char rhs);
    friend BasicSequence operator+(char lhs, BasicSequence&& rhs);
    friend BasicSequence operator+(BasicSequence&& lhs, char rhs);

    friend bool operator==(const BasicSequence& lhs, const BasicSequence& rhs);
    friend bool operator==(const std::string& lhs, const BasicSequence& rhs);
    friend bool operator==(const BasicSequence& lhs, const std::string& rhs);
    friend bool operator==(const char* lhs, const BasicSequence& rhs);
    friend bool operator==(const BasicSequence& lhs, const char* rhs);

    friend bool operator!=(const BasicSequence& lhs, const BasicSequence& rhs);
    friend bool operator!=(const std::string& lhs, const BasicSequence& rhs);
    friend bool operator!=(const BasicSequence& lhs, const std::string& rhs);
    friend bool operator!=(const char* lhs, const BasicSequence& rhs);
    friend bool operator!=(const BasicSequence& lhs, const char* rhs);

    friend bool operator<(const BasicSequence& lhs, const BasicSequence& rhs);
    friend bool operator<(const std::string& lhs, const BasicSequence& rhs);
    friend bool operator<(const BasicSequence& lhs, const std::string& rhs);
    friend bool operator<(const char* lhs, const BasicSequence& rhs);
    friend bool operator<(const BasicSequence& lhs, const char* rhs);

    friend bool operator<=(const BasicSequence& lhs, const BasicSequence& rhs);
    friend bool operator<=(const std::string& lhs, const BasicSequence& rhs);
    friend bool operator<=(const BasicSequence& lhs, const std::string& rhs);
    friend bool operator<=(const char* lhs, const BasicSequence& rhs);
    friend bool operator<=(const BasicSequence& lhs, const char* rhs);

    friend bool operator>(const BasicSequence& lhs, const BasicSequence& rhs);
    friend bool operator>(const std::string& lhs, const BasicSequence& rhs);
    friend bool operator>(const BasicSequence& lhs, const std::string& rhs);
    friend bool operator>(const char* lhs, const BasicSequence& rhs);
    friend bool operator>(const BasicSequence& lhs, const char* rhs);

    friend bool operator>=(const BasicSequence& lhs, const BasicSequence& rhs);
    friend bool operator>=(const std::string& lhs, const BasicSequence& rhs);
    friend bool operator>=(const BasicSequence& lhs, const std::string& rhs);
    friend bool operator>=(const char* lhs, const BasicSequence& rhs);
    friend bool operator>=(const BasicSequence& lhs, const char* rhs);

    friend void swap(BasicSequence& x, BasicSequence& y);

    friend std::istream& operator>>(std::istream& is, BasicSequence& rhs);
    friend std::ostream& operator<<(std::ostream& os, const BasicSequence& rhs);

    friend std::istream& getline(std::istream& is, BasicSequence& seq, char delim);
    friend std::istream& getline(std::istream&& is, BasicSequence& seq, char delim);
    friend std::istream& getline(std::istream& is, BasicSequence& seq);
    friend std::istream& getline(std::istream&& is, BasicSequence& seq);

    static void validate(const std::string& bases);
    static void validate(const char* bases);
    static void validate(std::initializer_list<char> bases);

protected:

    // No validation
    void set(std::string bases);
    void set(const char* bases);

    std::string s;

};

inline BasicSequence::Base& BasicSequence::Base::operator=(char base) { validate(base); b = base; return *this; }
inline BasicSequence::Base& BasicSequence::Base::operator=(Base base) { b = base.b; return *this; }
        
inline BasicSequence::Base::operator char() const { return b; }

inline void BasicSequence::Base::complement() { b = COMPLEMENTS[(unsigned char)b]; }
inline char BasicSequence::Base::operator~() const { return COMPLEMENTS[(unsigned char)b]; }

inline void BasicSequence::Base::validate(char base) { assert(COMPLEMENTS[(unsigned char)base]); }

inline BasicSequence::Base::Base(char& base): b(base) { validate(base); }
inline BasicSequence::Base::Base(Base& base): b(base.b) {}

const inline char BasicSequence::Base::COMPLEMENTS[256] = {
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//       !    "    #    $    %    &    '    (    )    *    +    ,    -    .    /
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , '.',  0 ,

//  0    1    2    3    4    5    6    7    8    9    :    ;    <    =    >    ?
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//  @    A    B    C    D    E    F    G    H    I    J    K    L    M    N    O
    0 , 'T', 'V', 'G', 'H',  0 ,  0 , 'C', 'D',  0 ,  0 , 'M',  0 , 'K', 'N',  0 , 

//  P    Q    R    S    T    U    V    W    X    Y    Z    [    \    ]    ^    _
    0 ,  0 , 'Y', 'S', 'A',  0 , 'B', 'W',  0 , 'R',  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//  `    a    b    c    d    e    f    g    h    i    j    k    l    m    n    o
    0 , 't', 'v', 'g', 'h',  0 ,  0 , 'c', 'd',  0 ,  0 , 'm',  0 , 'k', 'n',  0 ,

//  p    q    r    s    t    u    v    w    x    y    z    {    |    }    ~   DEL
    0 ,  0 , 'y', 's', 'a',  0 , 'b', 'w',  0 , 'r',  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0
};

inline BasicSequence::BasicSequence() = default;
inline BasicSequence::BasicSequence(const BasicSequence& seq): s(seq.s) {}
inline BasicSequence::BasicSequence(std::string bases): s(std::move(bases)) { validate(s); }
inline BasicSequence::BasicSequence(const BasicSequence& seq, size_t pos, size_t len): s(seq.s, pos, len) {}
inline BasicSequence::BasicSequence(const std::string& bases, size_t pos, size_t len): s(bases, pos, len) { validate(s); }
inline BasicSequence::BasicSequence(const char* bases): s(bases) { validate(s); }
inline BasicSequence::BasicSequence(const char* bases, size_t n): s(bases, n) { validate(s); }
inline BasicSequence::BasicSequence(size_t n, char base): s(n, base) { Base::validate(base); }
template <class InputIterator>
inline BasicSequence::BasicSequence(InputIterator first, InputIterator last): s(first, last) { validate(s); }
inline BasicSequence::BasicSequence(std::initializer_list<char> bases): s(bases) { validate(s); }
inline BasicSequence::BasicSequence(BasicSequence&& seq): s(seq.s) {}

inline BasicSequence& BasicSequence::operator=(const BasicSequence& seq) { s = seq.s; return *this; }
inline BasicSequence& BasicSequence::operator=(std::string bases) { validate(bases); s = std::move(bases); return *this; }
inline BasicSequence& BasicSequence::operator=(const char* bases) { validate(bases); s = bases; return *this; }
inline BasicSequence& BasicSequence::operator=(std::initializer_list<char> bases) { validate(bases); s = bases; return *this; }
inline BasicSequence& BasicSequence::operator=(BasicSequence&& seq) { s = seq.s; return *this; }
inline BasicSequence& BasicSequence::operator=(std::string&& bases) { validate(bases); s = bases; return *this; }

inline BasicSequence::operator const std::string&() const noexcept { return s; }
inline BasicSequence::operator const char*() const noexcept { return s.c_str(); }

inline BasicSequence::iterator BasicSequence::begin() noexcept { return s.begin(); }
inline BasicSequence::const_iterator BasicSequence::begin() const noexcept { return s.begin(); }
inline BasicSequence::iterator BasicSequence::end() noexcept { return s.end(); }
inline BasicSequence::const_iterator BasicSequence::end() const noexcept { return s.end(); }
inline BasicSequence::reverse_iterator BasicSequence::rbegin() noexcept { return s.rbegin(); }
inline BasicSequence::const_reverse_iterator BasicSequence::rbegin() const noexcept { return s.rbegin(); }
inline BasicSequence::reverse_iterator BasicSequence::rend() noexcept { return s.rend(); }
inline BasicSequence::const_reverse_iterator BasicSequence::rend() const noexcept { return s.rend(); }
inline BasicSequence::const_iterator BasicSequence::cbegin() const noexcept { return s.cbegin(); }
inline BasicSequence::const_iterator BasicSequence::cend() const noexcept { return s.cend(); }
inline BasicSequence::const_reverse_iterator BasicSequence::crbegin() const noexcept { return s.crbegin(); }
inline BasicSequence::const_reverse_iterator BasicSequence::crend() const noexcept { return s.crend(); }

inline size_t BasicSequence::size() const noexcept { return s.size(); }
inline size_t BasicSequence::length() const noexcept { return s.length(); }
inline size_t BasicSequence::max_size() const noexcept { return s.max_size(); }
inline void BasicSequence::resize (size_t n) { s.resize(n); }
inline void BasicSequence::resize (size_t n, char base) { Base::validate(base); s.resize(n, base); }
inline size_t BasicSequence::capacity() const noexcept { return s.capacity(); }
inline void BasicSequence::reserve(size_t n) { s.reserve(n); }
inline void BasicSequence::clear() noexcept { s.clear(); }
inline bool BasicSequence::empty() const noexcept { return s.empty(); }
inline void BasicSequence::shrink_to_fit() { s.shrink_to_fit(); }

inline BasicSequence::Base BasicSequence::operator[](size_t pos) { return Base(s[pos]); }
inline const char& BasicSequence::operator[](size_t pos) const { return s[pos]; }
inline BasicSequence::Base BasicSequence::at(size_t pos) { return Base(s.at(pos)); }
inline const char& BasicSequence::at(size_t pos) const { return s.at(pos); }
inline BasicSequence::Base BasicSequence::back() { return Base(s.back()); }
inline const char& BasicSequence::back() const { return s.back(); }
inline BasicSequence::Base BasicSequence::front() { return Base(s.front()); }
inline const char& BasicSequence::front() const { return s.front(); }

inline BasicSequence& BasicSequence::operator+=(const BasicSequence& rhs) { s += rhs.s; return *this; }
inline BasicSequence& BasicSequence::operator+=(const std::string& rhs) { validate(rhs); s += rhs; return *this; }
inline BasicSequence& BasicSequence::operator+=(const char* rhs) { validate(rhs); s += rhs; return *this; }
inline BasicSequence& BasicSequence::operator+=(char rhs) { Base::validate(rhs); s += rhs; return *this; }
inline BasicSequence& BasicSequence::operator+=(std::initializer_list<char> rhs) { validate(rhs); s += rhs; return *this; }

inline BasicSequence& BasicSequence::append(const BasicSequence& seq) { s.append(seq.s); return *this; }
inline BasicSequence& BasicSequence::append(const BasicSequence& seq, size_t subpos, size_t sublen) { s.append(seq.s, subpos, sublen); return *this; }
inline BasicSequence& BasicSequence::append(const std::string& bases) { validate(bases); s.append(bases); return *this; }
inline BasicSequence& BasicSequence::append(const std::string& bases, size_t subpos, size_t sublen) { validate(bases); s.append(bases, subpos, sublen); return *this; }
inline BasicSequence& BasicSequence::append(const char* bases) { validate(bases); s.append(bases); return *this; }
inline BasicSequence& BasicSequence::append(const char* bases, size_t n) { validate(bases); s.append(bases, n); return *this; }
inline BasicSequence& BasicSequence::append(size_t n, char base) { Base::validate(base); s.append(n, base); return *this; }
template <class InputIterator>
inline BasicSequence& BasicSequence::append(InputIterator first, InputIterator last) { s.append(first, last); validate(s); return *this; }
inline BasicSequence& BasicSequence::append(std::initializer_list<char> bases) { validate(bases); s.append(bases); return *this; }
inline void BasicSequence::push_back(char base) { Base::validate(base); s.push_back(base); }
inline BasicSequence& BasicSequence::assign(const BasicSequence& seq) { s.assign(seq.s); return *this; }
inline BasicSequence& BasicSequence::assign(const BasicSequence& seq, size_t subpos, size_t sublen) { s.assign(seq.s, subpos, sublen); return *this; }
inline BasicSequence& BasicSequence::assign(const std::string& bases) { validate(bases); s.assign(bases); return *this; }
inline BasicSequence& BasicSequence::assign(const std::string& bases, size_t subpos, size_t sublen) { validate(bases); s.assign(bases, subpos, npos); return *this; }
inline BasicSequence& BasicSequence::assign(const char* bases) { validate(bases); s.assign(bases); return *this; }
inline BasicSequence& BasicSequence::assign(const char* bases, size_t n) { validate(bases); s.assign(bases, n); return *this; }
inline BasicSequence& BasicSequence::assign(size_t n, char base) { Base::validate(base); s.assign(n, base); return *this; }
template <class InputIterator>
inline BasicSequence& BasicSequence::assign(InputIterator first, InputIterator last) { s.assign(first, last); validate(s); return *this; }
inline BasicSequence& BasicSequence::assign(std::initializer_list<char> bases) { validate(bases); s.assign(bases); return *this; }
inline BasicSequence& BasicSequence::assign(BasicSequence&& seq) noexcept { s.assign(seq.s); return *this; }
inline BasicSequence& BasicSequence::assign(std::string&& bases) noexcept { validate(bases); s.assign(bases); return *this; }

inline BasicSequence& BasicSequence::insert(size_t pos, const BasicSequence& seq) { s.insert(pos, seq.s); return *this; }
inline BasicSequence& BasicSequence::insert(size_t pos, const BasicSequence& seq, size_t subpos, size_t sublen) { s.insert(pos, seq.s, subpos, sublen); return *this; }
inline BasicSequence& BasicSequence::insert(size_t pos, const std::string& bases) { validate(bases); s.insert(pos, bases); return *this; }
inline BasicSequence& BasicSequence::insert(size_t pos, const std::string& bases, size_t subpos, size_t sublen) { validate(bases); s.insert(pos, bases, subpos, sublen); return *this; }
inline BasicSequence& BasicSequence::insert(size_t pos, const char* bases) { validate(bases); s.insert(pos, bases); return *this; }
inline BasicSequence& BasicSequence::insert(size_t pos, const char* bases, size_t n) { validate(bases); s.insert(pos,bases, n); return *this; }
inline BasicSequence& BasicSequence::insert(size_t pos, size_t n, char base) { Base::validate(base); s.insert(pos, n, base); return *this; }
inline BasicSequence::iterator BasicSequence::insert(const_iterator p, size_t n, char base) { Base::validate(base); return s.insert(p, n, base); }
inline BasicSequence::iterator BasicSequence::insert(const_iterator p, char base) { Base::validate(base); return s.insert(p, base); }
template <class InputIterator>
inline BasicSequence::iterator BasicSequence::insert(iterator p, InputIterator first, InputIterator last) { auto ret = s.insert(p, first, last); validate(s); return ret; }
inline BasicSequence& BasicSequence::insert(const_iterator p, std::initializer_list<char> bases) { validate(bases); s.insert(p, bases); return *this; }

inline BasicSequence& BasicSequence::erase (size_t pos, size_t len) { s.erase(pos, len); return *this; }
inline BasicSequence::iterator BasicSequence::erase (iterator p) { return s.erase(p); }
inline BasicSequence::iterator BasicSequence::erase (iterator first, iterator last) { return s.erase(first, last); }

inline BasicSequence& BasicSequence::replace(size_t pos, size_t len, const BasicSequence& seq) { s.replace(pos, len, seq.s); return *this; }
inline BasicSequence& BasicSequence::replace(size_t pos, size_t len, const std::string& bases) { validate(bases); s.replace(pos, len, bases); return *this; }
inline BasicSequence& BasicSequence::replace(const_iterator i1, const_iterator i2, const BasicSequence& seq) { s.replace(i1, i2, seq.s); return *this; }
inline BasicSequence& BasicSequence::replace(const_iterator i1, const_iterator i2, const std::string& bases) { validate(bases); s.replace(i1, i2, bases); return *this; }
inline BasicSequence& BasicSequence::replace(size_t pos, size_t len, const BasicSequence& seq, size_t subpos, size_t sublen) { s.replace(pos, len, seq.s, subpos, sublen); return *this; }
inline BasicSequence& BasicSequence::replace(size_t pos, size_t len, const std::string& bases, size_t subpos, size_t sublen) { validate(bases); s.replace(pos, len, bases, subpos, sublen); return *this; }
inline BasicSequence& BasicSequence::replace(size_t pos, size_t len, const char* bases) { validate(bases); s.replace(pos, len, bases); return *this; }
inline BasicSequence& BasicSequence::replace(const_iterator i1, const_iterator i2, const char* bases) { validate(bases); s.replace(i1, i2, bases); return *this; }
inline BasicSequence& BasicSequence::replace(size_t pos, size_t len, const char* bases, size_t n) { validate(bases); s.replace(pos, len, bases, n); return *this; }
inline BasicSequence& BasicSequence::replace(const_iterator i1, const_iterator i2, const char* bases, size_t n) { validate(bases); s.replace(i1, i2, bases, n); return *this; }
inline BasicSequence& BasicSequence::replace(size_t pos, size_t len, size_t n, char base) { Base::validate(base); s.replace(pos, len, n, base); return *this; }
inline BasicSequence& BasicSequence::replace(const_iterator i1, const_iterator i2, size_t n, char base) { Base::validate(base); s.replace(i1, i2, n, base); return *this; }
template <class InputIterator>
inline BasicSequence& BasicSequence::replace(const_iterator i1, const_iterator i2, InputIterator first, InputIterator last) { s.replace(i1, i2, first, last); validate(s); return *this; }
inline BasicSequence& BasicSequence::replace(const_iterator i1, const_iterator i2, std::initializer_list<char> bases) { validate(bases); s.replace(i1, i2, bases); return *this; }

inline void BasicSequence::swap(BasicSequence& seq) { s.swap(seq.s); }
inline void BasicSequence::swap(std::string& bases) { validate(bases); s.swap(bases); }

inline void BasicSequence::pop_back() { s.pop_back(); }

inline const char* BasicSequence::c_str() const noexcept { return s.c_str(); }
inline const char* BasicSequence::data() const noexcept { return s.data(); }

inline BasicSequence::allocator_type BasicSequence::get_allocator() const noexcept { return s.get_allocator(); }

inline size_t BasicSequence::copy(char* bases, size_t len, size_t pos) const { return s.copy(bases, len, pos); }

inline size_t BasicSequence::find(const BasicSequence& seq, size_t pos) const noexcept { return s.find(seq.s, pos); }
inline size_t BasicSequence::find(const std::string& bases, size_t pos) const noexcept { return s.find(bases, pos); }
inline size_t BasicSequence::find(const char* bases, size_t pos) const { return s.find(bases, pos); }
inline size_t BasicSequence::find(const char* bases, size_t pos, size_type n) const { return s.find(bases, pos, n); }
inline size_t BasicSequence::find(char base, size_t pos) const noexcept { return s.find(base, pos); }

inline size_t BasicSequence::rfind(const BasicSequence& seq, size_t pos) const noexcept { return s.rfind(seq.s, pos); }
inline size_t BasicSequence::rfind(const std::string& bases, size_t pos) const noexcept { return s.rfind(bases, pos); }
inline size_t BasicSequence::rfind(const char* bases, size_t pos) const { return s.rfind(bases, pos); }
inline size_t BasicSequence::rfind(const char* bases, size_t pos, size_t n) const { return s.rfind(bases, pos, n); }
inline size_t BasicSequence::rfind(char base, size_t pos) const noexcept { return s.rfind(base, pos); }

inline size_t BasicSequence::find_first_of(const BasicSequence& seq, size_t pos) const noexcept { return s.find_first_of(seq.s, pos); }
inline size_t BasicSequence::find_first_of(const std::string& bases, size_t pos) const noexcept { return s.find_first_of(bases, pos); }
inline size_t BasicSequence::find_first_of(const char* bases, size_t pos) const { return s.find_first_of(bases, pos); }
inline size_t BasicSequence::find_first_of(const char* bases, size_t pos, size_t n) const { return s.find_first_of(bases, pos, n); }
inline size_t BasicSequence::find_first_of(char base, size_t pos) const noexcept { return s.find_first_of(base, pos); }

inline size_t BasicSequence::find_last_of(const BasicSequence& seq, size_t pos) const noexcept { return s.find_last_of(seq.s, pos); }
inline size_t BasicSequence::find_last_of(const std::string& bases, size_t pos) const noexcept { return s.find_last_of(bases, pos); }
inline size_t BasicSequence::find_last_of(const char* bases, size_t pos) const { return s.find_last_of(bases, pos); }
inline size_t BasicSequence::find_last_of(const char* bases, size_t pos, size_t n) const { return s.find_last_of(bases, pos, n); }
inline size_t BasicSequence::find_last_of(char base, size_t pos) const noexcept { return s.find_last_of(base, pos); }

inline size_t BasicSequence::find_first_not_of(const BasicSequence& seq, size_t pos) const noexcept { return s.find_first_not_of(seq.s, pos); }
inline size_t BasicSequence::find_first_not_of(const std::string& bases, size_t pos) const noexcept { return s.find_first_not_of(bases, pos); }
inline size_t BasicSequence::find_first_not_of(const char* bases, size_t pos) const { return s.find_first_not_of(bases, pos); }
inline size_t BasicSequence::find_first_not_of(const char* bases, size_t pos, size_t n) const { return s.find_first_not_of(bases, pos, n); }
inline size_t BasicSequence::find_first_not_of(char base, size_t pos) const noexcept { return s.find_first_not_of(base, pos); }

inline size_t BasicSequence::find_last_not_of(const BasicSequence& seq, size_t pos) const noexcept { return s.find_last_not_of(seq.s, pos); }
inline size_t BasicSequence::find_last_not_of(const std::string& bases, size_t pos) const noexcept { return s.find_last_not_of(bases, pos); }
inline size_t BasicSequence::find_last_not_of(const char* bases, size_t pos) const { return s.find_last_not_of(bases, pos); }
inline size_t BasicSequence::find_last_not_of(const char* bases, size_t pos, size_t n) const { return s.find_last_not_of(bases, pos, n); }
inline size_t BasicSequence::find_last_not_of(char base, size_t pos) const noexcept { return s.find_last_not_of(base, pos); }

inline BasicSequence BasicSequence::substr(size_t pos, size_t len) const { BasicSequence seq; seq.set(s.substr(pos, len)); return seq; }

inline int BasicSequence::compare(const BasicSequence& seq) const noexcept { return s.compare(seq.s); }
inline int BasicSequence::compare(const std::string& bases) const noexcept { return s.compare(bases); }
inline int BasicSequence::compare(size_t pos, size_t len, const BasicSequence& seq) const { return s.compare(pos, len, seq.s); }
inline int BasicSequence::compare(size_t pos, size_t len, const std::string& bases) const { return s.compare(pos, len, bases); }
inline int BasicSequence::compare(size_t pos, size_t len, const BasicSequence& seq, size_t subpos, size_t sublen) const { return s.compare(pos, len, seq.s, subpos, sublen); }
inline int BasicSequence::compare(size_t pos, size_t len, const std::string& bases, size_t subpos, size_t sublen) const { return s.compare(pos, len, bases, subpos, sublen); }
inline int BasicSequence::compare(const char* bases) const { return s.compare(bases); }
inline int BasicSequence::compare(size_t pos, size_t len, const char* bases) const { return s.compare(pos, len, bases); }
inline int BasicSequence::compare(size_t pos, size_t len, const char* bases, size_t n) const { return s.compare(pos, len, bases, n); }

inline BasicSequence operator+(const BasicSequence& lhs, const BasicSequence& rhs) { return lhs.s + rhs.s; }
inline BasicSequence operator+(BasicSequence&& lhs, BasicSequence&& rhs) { return lhs.s + rhs.s; }
inline BasicSequence operator+(BasicSequence&& lhs, const BasicSequence& rhs) { return lhs.s + rhs.s; }
inline BasicSequence operator+(const BasicSequence& lhs, BasicSequence&& rhs) { return lhs.s + rhs.s; }
inline BasicSequence operator+(const std::string& lhs, const BasicSequence& rhs) { return lhs + rhs.s; }
inline BasicSequence operator+(const BasicSequence& lhs, const std::string& rhs) { return lhs.s + rhs; }
inline BasicSequence operator+(const std::string& lhs, BasicSequence&& rhs) { return lhs + rhs.s; }
inline BasicSequence operator+(const BasicSequence& lhs, std::string&& rhs) { return lhs.s + rhs; }
inline BasicSequence operator+(std::string&& lhs, const BasicSequence& rhs) { return lhs + rhs.s; }
inline BasicSequence operator+(BasicSequence&& lhs, const std::string& rhs) { return lhs.s + rhs; }
inline BasicSequence operator+(BasicSequence&& lhs, std::string&& rhs) { return lhs.s + rhs; }
inline BasicSequence operator+(std::string&& lhs, BasicSequence&& rhs) { return lhs + rhs.s; }
inline BasicSequence operator+(const char* lhs, const BasicSequence& rhs) { return lhs + rhs.s; }
inline BasicSequence operator+(const BasicSequence& lhs, const char* rhs) { return lhs.s + rhs; }
inline BasicSequence operator+(const char* lhs, BasicSequence&& rhs) { return lhs + rhs.s; }
inline BasicSequence operator+(BasicSequence&& lhs, const char* rhs) { return lhs.s + rhs; }
inline BasicSequence operator+(char lhs, const BasicSequence& rhs) { return lhs + rhs.s; }
inline BasicSequence operator+(const BasicSequence& lhs, char rhs) { return lhs.s + rhs; }
inline BasicSequence operator+(char lhs, BasicSequence&& rhs) { return lhs + rhs.s; }
inline BasicSequence operator+(BasicSequence&& lhs, char rhs) { return lhs.s + rhs; }

inline bool operator==(const BasicSequence& lhs, const BasicSequence& rhs) { return lhs.s == rhs.s; }
inline bool operator==(const std::string& lhs, const BasicSequence& rhs) { return lhs == rhs.s; }
inline bool operator==(const BasicSequence& lhs, const std::string& rhs) { return lhs.s == rhs; }
inline bool operator==(const char* lhs, const BasicSequence& rhs) { return lhs == rhs.s; }
inline bool operator==(const BasicSequence& lhs, const char* rhs) { return lhs.s == rhs; }

inline bool operator!=(const BasicSequence& lhs, const BasicSequence& rhs) { return lhs.s != rhs.s; }
inline bool operator!=(const std::string& lhs, const BasicSequence& rhs) { return lhs != rhs.s; }
inline bool operator!=(const BasicSequence& lhs, const std::string& rhs) { return lhs.s != rhs; }
inline bool operator!=(const char* lhs, const BasicSequence& rhs) { return lhs != rhs.s; }
inline bool operator!=(const BasicSequence& lhs, const char* rhs) { return lhs.s != rhs; }

inline bool operator<(const BasicSequence& lhs, const BasicSequence& rhs) { return lhs.s < rhs.s; }
inline bool operator<(const std::string& lhs, const BasicSequence& rhs) { return lhs < rhs.s; }
inline bool operator<(const BasicSequence& lhs, const std::string& rhs) { return lhs.s < rhs; }
inline bool operator<(const char* lhs, const BasicSequence& rhs) { return lhs < rhs.s; }
inline bool operator<(const BasicSequence& lhs, const char* rhs) { return lhs.s < rhs; }

inline bool operator<=(const BasicSequence& lhs, const BasicSequence& rhs) { return lhs.s <= rhs.s; }
inline bool operator<=(const std::string& lhs, const BasicSequence& rhs) { return lhs <= rhs.s; }
inline bool operator<=(const BasicSequence& lhs, const std::string& rhs) { return lhs.s <= rhs; }
inline bool operator<=(const char* lhs, const BasicSequence& rhs) { return lhs <= rhs.s; }
inline bool operator<=(const BasicSequence& lhs, const char* rhs) { return lhs.s <= rhs; }

inline bool operator>(const BasicSequence& lhs, const BasicSequence& rhs) { return lhs.s > rhs.s; }
inline bool operator>(const std::string& lhs, const BasicSequence& rhs) { return lhs > rhs.s; }
inline bool operator>(const BasicSequence& lhs, const std::string& rhs) { return lhs.s > rhs; }
inline bool operator>(const char* lhs, const BasicSequence& rhs) { return lhs > rhs.s; }
inline bool operator>(const BasicSequence& lhs, const char* rhs) { return lhs.s > rhs; }

inline bool operator>=(const BasicSequence& lhs, const BasicSequence& rhs) { return lhs.s >= rhs.s; }
inline bool operator>=(const std::string& lhs, const BasicSequence& rhs) { return lhs >= rhs.s; }
inline bool operator>=(const BasicSequence& lhs, const std::string& rhs) { return lhs.s >= rhs; }
inline bool operator>=(const char* lhs, const BasicSequence& rhs) { return lhs >= rhs.s; }
inline bool operator>=(const BasicSequence& lhs, const char* rhs) { return lhs.s >= rhs; }

inline void swap(BasicSequence& x, BasicSequence& y) { swap(x.s, y.s); }

inline std::istream& operator>>(std::istream& is, BasicSequence& rhs) { is >> rhs.s; return is; }
inline std::ostream& operator<<(std::ostream& os, const BasicSequence& rhs) { os << rhs.s; return os; }

inline std::istream& getline(std::istream& is, BasicSequence& seq, char delim) { auto& ret = getline(is, seq.s, delim); BasicSequence::validate(seq.s); return ret; }
inline std::istream& getline(std::istream&& is, BasicSequence& seq, char delim) { auto& ret = getline(is, seq.s, delim); BasicSequence::validate(seq.s); return ret; }
inline std::istream& getline(std::istream& is, BasicSequence& seq) { auto& ret = getline(is, seq.s); BasicSequence::validate(seq.s); return ret; }
inline std::istream& getline(std::istream&& is, BasicSequence& seq) { auto& ret = getline(is, seq.s); BasicSequence::validate(seq.s); return ret; }

inline void BasicSequence::validate(const std::string& bases) {
    std::for_each(bases.begin(), bases.end(),
        [] (char base) { Base::validate(base); }
    );
}

inline void BasicSequence::validate(const char* bases) {
    for (unsigned i = 0; bases[i] != '\0'; ++i) {
        Base::validate(bases[i]);  
    }
}

inline void BasicSequence::validate(std::initializer_list<char> bases) {
    std::for_each(bases.begin(), bases.end(),
        [] (char base) { Base::validate(base); }
    );
}

inline void BasicSequence::set(std::string bases) { s = std::move(bases); }
inline void BasicSequence::set(const char* bases) { s = std::string(bases); }

} // namespace btl

#endif
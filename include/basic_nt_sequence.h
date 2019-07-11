#ifndef BASIC_NT_SEQUENCE_H
#define BASIC_NT_SEQUENCE_H

#include <string>
#include <algorithm>
#include <cassert>
#include <iostream>

namespace btl {

class BasicNtSequence {

public:

    class Base {

    public:

        Base& operator=(char base);
        Base& operator=(Base base);

        operator char() const;

        void complement();
        char operator~() const;

        static inline void validate(char base);

        friend class BasicNtSequence;

        Base(char& base);
        Base(Base& base);

        char& b;
        static const char COMPLEMENTS[256];
    };

    static inline const size_t& npos = std::string::npos;

    BasicNtSequence();
    BasicNtSequence(const BasicNtSequence& seq);
    BasicNtSequence(std::string bases);
    BasicNtSequence(const BasicNtSequence& seq, size_t pos, size_t len = npos);
    BasicNtSequence(const std::string& bases, size_t pos, size_t len = npos);
    BasicNtSequence(const char* bases);
    BasicNtSequence(const char* bases, size_t n);
    BasicNtSequence(size_t n, char base);
    template <class InputIterator>
    BasicNtSequence(InputIterator first, InputIterator last);
    BasicNtSequence(std::initializer_list<char> bases);
    BasicNtSequence(BasicNtSequence&& seq);

    BasicNtSequence& operator=(const BasicNtSequence& seq);
    BasicNtSequence& operator=(std::string bases);
    BasicNtSequence& operator=(const char* bases);
    BasicNtSequence& operator=(std::initializer_list<char> bases);
    BasicNtSequence& operator=(BasicNtSequence&& seq);
    BasicNtSequence& operator=(std::string&& bases);

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

    BasicNtSequence& operator+=(const BasicNtSequence& rhs);
    BasicNtSequence& operator+=(const std::string& rhs);
    BasicNtSequence& operator+=(const char* rhs);
    BasicNtSequence& operator+=(char rhs);
    BasicNtSequence& operator+=(std::initializer_list<char> rhs);

    BasicNtSequence& append(const BasicNtSequence& seq);
    BasicNtSequence& append(const BasicNtSequence& seq, size_t subpos, size_t sublen);
    BasicNtSequence& append(const std::string& bases);
    BasicNtSequence& append(const std::string& bases, size_t subpos, size_t sublen);
    BasicNtSequence& append(const char* bases);
    BasicNtSequence& append(const char* bases, size_t n);
    BasicNtSequence& append(size_t n, char base);
    template <class InputIterator>
    BasicNtSequence& append(InputIterator first, InputIterator last);
    BasicNtSequence& append(std::initializer_list<char> bases);

    void push_back(char base);

    BasicNtSequence& assign(const BasicNtSequence& seq);
    BasicNtSequence& assign(const BasicNtSequence& seq, size_t subpos, size_t sublen = npos);
    BasicNtSequence& assign(const std::string& bases);
    BasicNtSequence& assign(const std::string& bases, size_t subpos, size_t sublen = npos);
    BasicNtSequence& assign(const char* bases);
    BasicNtSequence& assign(const char* bases, size_t n);
    BasicNtSequence& assign(size_t n, char base);
    template <class InputIterator>
    BasicNtSequence& assign(InputIterator first, InputIterator last);
    BasicNtSequence& assign(std::initializer_list<char> bases);
    BasicNtSequence& assign(BasicNtSequence&& seq) noexcept;
    BasicNtSequence& assign(std::string&& bases) noexcept;

    BasicNtSequence& insert(size_t pos, const BasicNtSequence& seq);
    BasicNtSequence& insert(size_t pos, const BasicNtSequence& seq, size_t subpos, size_t sublen = npos);
    BasicNtSequence& insert(size_t pos, const std::string& bases);
    BasicNtSequence& insert(size_t pos, const std::string& bases, size_t subpos, size_t sublen = npos);
    BasicNtSequence& insert(size_t pos, const char* bases);
    BasicNtSequence& insert(size_t pos, const char* bases, size_t n);
    BasicNtSequence& insert(size_t pos, size_t n, char base);
    iterator insert(const_iterator p, size_t n, char base);
    iterator insert(const_iterator p, char base);
    template <class InputIterator>
    iterator insert(iterator p, InputIterator first, InputIterator last);
    BasicNtSequence& insert(const_iterator p, std::initializer_list<char> bases);

    BasicNtSequence& erase(size_t pos = 0, size_t len = npos);
    iterator erase(iterator p);
    iterator erase(iterator first, iterator last);

    BasicNtSequence& replace(size_t pos, size_t len, const BasicNtSequence& seq);
    BasicNtSequence& replace(size_t pos, size_t len, const std::string& bases);
    BasicNtSequence& replace(const_iterator i1, const_iterator i2, const BasicNtSequence& seq);
    BasicNtSequence& replace(const_iterator i1, const_iterator i2, const std::string& bases);
    BasicNtSequence& replace(size_t pos, size_t len, const BasicNtSequence& seq, size_t subpos, size_t sublen = npos);
    BasicNtSequence& replace(size_t pos, size_t len, const std::string& bases, size_t subpos, size_t sublen = npos);
    BasicNtSequence& replace(size_t pos, size_t len, const char* bases);
    BasicNtSequence& replace(const_iterator i1, const_iterator i2, const char* bases);
    BasicNtSequence& replace(size_t pos, size_t len, const char* bases, size_t n);
    BasicNtSequence& replace(const_iterator i1, const_iterator i2, const char* bases, size_t n);
    BasicNtSequence& replace(size_t pos, size_t len, size_t n, char base);
    BasicNtSequence& replace(const_iterator i1, const_iterator i2, size_t n, char base);
    template <class InputIterator>
    BasicNtSequence& replace(const_iterator i1, const_iterator i2, InputIterator first, InputIterator last);
    BasicNtSequence& replace(const_iterator i1, const_iterator i2, std::initializer_list<char> bases);

    void swap(BasicNtSequence& seq);
    void swap(std::string& bases);

    void pop_back();

    const char* c_str() const noexcept;
    const char* data() const noexcept;

    using allocator_type = std::string::allocator_type;

    allocator_type get_allocator() const noexcept;

    size_t copy(char* bases, size_t len, size_t pos = 0) const;

    using size_type = std::string::size_type;

    size_t find(const BasicNtSequence& seq, size_t pos = 0) const noexcept;
    size_t find(const std::string& bases, size_t pos = 0) const noexcept;
    size_t find(const char* bases, size_t pos = 0) const;
    size_t find(const char* bases, size_t pos, size_type n) const;
    size_t find(char base, size_t pos = 0) const noexcept;

    size_t rfind(const BasicNtSequence& seq, size_t pos = npos) const noexcept;
    size_t rfind(const std::string& bases, size_t pos = npos) const noexcept;
    size_t rfind(const char* bases, size_t pos = npos) const;
    size_t rfind(const char* bases, size_t pos, size_t n) const;
    size_t rfind(char base, size_t pos = npos) const noexcept;

    size_t find_first_of(const BasicNtSequence& seq, size_t pos = 0) const noexcept;
    size_t find_first_of(const std::string& bases, size_t pos = 0) const noexcept;
    size_t find_first_of(const char* bases, size_t pos = 0) const;
    size_t find_first_of(const char* bases, size_t pos, size_t n) const;
    size_t find_first_of(char base, size_t pos = 0) const noexcept;

    size_t find_last_of(const BasicNtSequence& seq, size_t pos = npos) const noexcept;
    size_t find_last_of(const std::string& bases, size_t pos = npos) const noexcept;
    size_t find_last_of(const char* bases, size_t pos = npos) const;
    size_t find_last_of(const char* bases, size_t pos, size_t n) const;
    size_t find_last_of(char base, size_t pos = npos) const noexcept;

    size_t find_first_not_of(const BasicNtSequence& seq, size_t pos = 0) const noexcept;
    size_t find_first_not_of(const std::string& bases, size_t pos = 0) const noexcept;
    size_t find_first_not_of(const char* bases, size_t pos = 0) const;
    size_t find_first_not_of(const char* bases, size_t pos, size_t n) const;
    size_t find_first_not_of(char base, size_t pos = 0) const noexcept;

    size_t find_last_not_of(const BasicNtSequence& seq, size_t pos = npos) const noexcept;
    size_t find_last_not_of(const std::string& bases, size_t pos = npos) const noexcept;
    size_t find_last_not_of(const char* bases, size_t pos = npos) const;
    size_t find_last_not_of(const char* bases, size_t pos, size_t n) const;
    size_t find_last_not_of(char base, size_t pos = npos) const noexcept;

    BasicNtSequence substr(size_t pos = 0, size_t len = npos) const;

    int compare(const BasicNtSequence& seq) const noexcept;
    int compare(const std::string& bases) const noexcept;
    int compare(size_t pos, size_t len, const BasicNtSequence& seq) const;
    int compare(size_t pos, size_t len, const std::string& bases) const;
    int compare(size_t pos, size_t len, const BasicNtSequence& seq, size_t subpos, size_t sublen = npos) const;
    int compare(size_t pos, size_t len, const std::string& bases, size_t subpos, size_t sublen = npos) const;
    int compare(const char* bases) const;
    int compare(size_t pos, size_t len, const char* bases) const;
    int compare(size_t pos, size_t len, const char* bases, size_t n) const;

    friend BasicNtSequence operator+(const BasicNtSequence& lhs, const BasicNtSequence& rhs);
    friend BasicNtSequence operator+(BasicNtSequence&& lhs, BasicNtSequence&& rhs);
    friend BasicNtSequence operator+(BasicNtSequence&& lhs, const BasicNtSequence& rhs);
    friend BasicNtSequence operator+(const BasicNtSequence& lhs, BasicNtSequence&& rhs);
    friend BasicNtSequence operator+(const std::string& lhs, const BasicNtSequence& rhs);
    friend BasicNtSequence operator+(const BasicNtSequence& lhs, const std::string& rhs);
    friend BasicNtSequence operator+(const std::string& lhs, BasicNtSequence&& rhs);
    friend BasicNtSequence operator+(const BasicNtSequence& lhs, std::string&& rhs);
    friend BasicNtSequence operator+(std::string&& lhs, const BasicNtSequence& rhs);
    friend BasicNtSequence operator+(BasicNtSequence&& lhs, const std::string& rhs);
    friend BasicNtSequence operator+(BasicNtSequence&& lhs, std::string&& rhs);
    friend BasicNtSequence operator+(std::string&& lhs, BasicNtSequence&& rhs);
    friend BasicNtSequence operator+(const char* lhs, const BasicNtSequence& rhs);
    friend BasicNtSequence operator+(const BasicNtSequence& lhs, const char* rhs);
    friend BasicNtSequence operator+(const char* lhs, BasicNtSequence&& rhs);
    friend BasicNtSequence operator+(BasicNtSequence&& lhs, const char* rhs);
    friend BasicNtSequence operator+(char lhs, const BasicNtSequence& rhs);
    friend BasicNtSequence operator+(const BasicNtSequence& lhs, char rhs);
    friend BasicNtSequence operator+(char lhs, BasicNtSequence&& rhs);
    friend BasicNtSequence operator+(BasicNtSequence&& lhs, char rhs);

    friend bool operator==(const BasicNtSequence& lhs, const BasicNtSequence& rhs);
    friend bool operator==(const std::string& lhs, const BasicNtSequence& rhs);
    friend bool operator==(const BasicNtSequence& lhs, const std::string& rhs);
    friend bool operator==(const char* lhs, const BasicNtSequence& rhs);
    friend bool operator==(const BasicNtSequence& lhs, const char* rhs);

    friend bool operator!=(const BasicNtSequence& lhs, const BasicNtSequence& rhs);
    friend bool operator!=(const std::string& lhs, const BasicNtSequence& rhs);
    friend bool operator!=(const BasicNtSequence& lhs, const std::string& rhs);
    friend bool operator!=(const char* lhs, const BasicNtSequence& rhs);
    friend bool operator!=(const BasicNtSequence& lhs, const char* rhs);

    friend bool operator<(const BasicNtSequence& lhs, const BasicNtSequence& rhs);
    friend bool operator<(const std::string& lhs, const BasicNtSequence& rhs);
    friend bool operator<(const BasicNtSequence& lhs, const std::string& rhs);
    friend bool operator<(const char* lhs, const BasicNtSequence& rhs);
    friend bool operator<(const BasicNtSequence& lhs, const char* rhs);

    friend bool operator<=(const BasicNtSequence& lhs, const BasicNtSequence& rhs);
    friend bool operator<=(const std::string& lhs, const BasicNtSequence& rhs);
    friend bool operator<=(const BasicNtSequence& lhs, const std::string& rhs);
    friend bool operator<=(const char* lhs, const BasicNtSequence& rhs);
    friend bool operator<=(const BasicNtSequence& lhs, const char* rhs);

    friend bool operator>(const BasicNtSequence& lhs, const BasicNtSequence& rhs);
    friend bool operator>(const std::string& lhs, const BasicNtSequence& rhs);
    friend bool operator>(const BasicNtSequence& lhs, const std::string& rhs);
    friend bool operator>(const char* lhs, const BasicNtSequence& rhs);
    friend bool operator>(const BasicNtSequence& lhs, const char* rhs);

    friend bool operator>=(const BasicNtSequence& lhs, const BasicNtSequence& rhs);
    friend bool operator>=(const std::string& lhs, const BasicNtSequence& rhs);
    friend bool operator>=(const BasicNtSequence& lhs, const std::string& rhs);
    friend bool operator>=(const char* lhs, const BasicNtSequence& rhs);
    friend bool operator>=(const BasicNtSequence& lhs, const char* rhs);

    friend void swap(BasicNtSequence& x, BasicNtSequence& y);

    friend std::istream& operator>>(std::istream& is, BasicNtSequence& rhs);
    friend std::ostream& operator<<(std::ostream& os, const BasicNtSequence& rhs);

    friend std::istream& getline(std::istream& is, BasicNtSequence& seq, char delim);
    friend std::istream& getline(std::istream&& is, BasicNtSequence& seq, char delim);
    friend std::istream& getline(std::istream& is, BasicNtSequence& seq);
    friend std::istream& getline(std::istream&& is, BasicNtSequence& seq);

    static void validate(const std::string& bases);
    static void validate(const char* bases);
    static void validate(std::initializer_list<char> bases);

protected:

    // No validation
    void set(std::string bases);
    void set(const char* bases);

    std::string s;

};

inline BasicNtSequence::Base& BasicNtSequence::Base::operator=(char base) { validate(base); b = base; return *this; }
inline BasicNtSequence::Base& BasicNtSequence::Base::operator=(Base base) { b = base.b; return *this; }
        
inline BasicNtSequence::Base::operator char() const { return b; }

inline void BasicNtSequence::Base::complement() { b = COMPLEMENTS[(unsigned char)b]; }
inline char BasicNtSequence::Base::operator~() const { return COMPLEMENTS[(unsigned char)b]; }

inline void BasicNtSequence::Base::validate(char base) { assert(COMPLEMENTS[(unsigned char)base]); }

inline BasicNtSequence::Base::Base(char& base): b(base) { validate(base); }
inline BasicNtSequence::Base::Base(Base& base): b(base.b) {}

const inline char BasicNtSequence::Base::COMPLEMENTS[256] = {
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//       !    "    #    $    %    &    '    (    )    *    +    ,    -    .    /
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , '-' , '.',  0 ,

//  0    1    2    3    4    5    6    7    8    9    :    ;    <    =    >    ?
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//  @    A    B    C    D    E    F    G    H    I    J    K    L    M    N    O
    0 , 'T', 'V', 'G', 'H',  0 ,  0 , 'C', 'D',  0 ,  0 , 'M',  0 , 'K', 'N',  0 , 

//  P    Q    R    S    T    U    V    W    X    Y    Z    [    \    ]    ^    _
    0 ,  0 , 'Y', 'S', 'A', 'U' , 'B', 'W',  0 , 'R',  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//  `    a    b    c    d    e    f    g    h    i    j    k    l    m    n    o
    0 , 't', 'v', 'g', 'h',  0 ,  0 , 'c', 'd',  0 ,  0 , 'm',  0 , 'k', 'n',  0 ,

//  p    q    r    s    t    u    v    w    x    y    z    {    |    }    ~   DEL
    0 ,  0 , 'y', 's', 'a', 'u' , 'b', 'w',  0 , 'r',  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0
};

inline BasicNtSequence::BasicNtSequence() = default;
inline BasicNtSequence::BasicNtSequence(const BasicNtSequence& seq): s(seq.s) {}
inline BasicNtSequence::BasicNtSequence(std::string bases): s(std::move(bases)) { validate(s); }
inline BasicNtSequence::BasicNtSequence(const BasicNtSequence& seq, size_t pos, size_t len): s(seq.s, pos, len) {}
inline BasicNtSequence::BasicNtSequence(const std::string& bases, size_t pos, size_t len): s(bases, pos, len) { validate(s); }
inline BasicNtSequence::BasicNtSequence(const char* bases): s(bases) { validate(s); }
inline BasicNtSequence::BasicNtSequence(const char* bases, size_t n): s(bases, n) { validate(s); }
inline BasicNtSequence::BasicNtSequence(size_t n, char base): s(n, base) { Base::validate(base); }
template <class InputIterator>
inline BasicNtSequence::BasicNtSequence(InputIterator first, InputIterator last): s(first, last) { validate(s); }
inline BasicNtSequence::BasicNtSequence(std::initializer_list<char> bases): s(bases) { validate(s); }
inline BasicNtSequence::BasicNtSequence(BasicNtSequence&& seq): s(seq.s) {}

inline BasicNtSequence& BasicNtSequence::operator=(const BasicNtSequence& seq) { s = seq.s; return *this; }
inline BasicNtSequence& BasicNtSequence::operator=(std::string bases) { validate(bases); s = std::move(bases); return *this; }
inline BasicNtSequence& BasicNtSequence::operator=(const char* bases) { validate(bases); s = bases; return *this; }
inline BasicNtSequence& BasicNtSequence::operator=(std::initializer_list<char> bases) { validate(bases); s = bases; return *this; }
inline BasicNtSequence& BasicNtSequence::operator=(BasicNtSequence&& seq) { s = seq.s; return *this; }
inline BasicNtSequence& BasicNtSequence::operator=(std::string&& bases) { validate(bases); s = bases; return *this; }

inline BasicNtSequence::operator const std::string&() const noexcept { return s; }
inline BasicNtSequence::operator const char*() const noexcept { return s.c_str(); }

inline BasicNtSequence::iterator BasicNtSequence::begin() noexcept { return s.begin(); }
inline BasicNtSequence::const_iterator BasicNtSequence::begin() const noexcept { return s.begin(); }
inline BasicNtSequence::iterator BasicNtSequence::end() noexcept { return s.end(); }
inline BasicNtSequence::const_iterator BasicNtSequence::end() const noexcept { return s.end(); }
inline BasicNtSequence::reverse_iterator BasicNtSequence::rbegin() noexcept { return s.rbegin(); }
inline BasicNtSequence::const_reverse_iterator BasicNtSequence::rbegin() const noexcept { return s.rbegin(); }
inline BasicNtSequence::reverse_iterator BasicNtSequence::rend() noexcept { return s.rend(); }
inline BasicNtSequence::const_reverse_iterator BasicNtSequence::rend() const noexcept { return s.rend(); }
inline BasicNtSequence::const_iterator BasicNtSequence::cbegin() const noexcept { return s.cbegin(); }
inline BasicNtSequence::const_iterator BasicNtSequence::cend() const noexcept { return s.cend(); }
inline BasicNtSequence::const_reverse_iterator BasicNtSequence::crbegin() const noexcept { return s.crbegin(); }
inline BasicNtSequence::const_reverse_iterator BasicNtSequence::crend() const noexcept { return s.crend(); }

inline size_t BasicNtSequence::size() const noexcept { return s.size(); }
inline size_t BasicNtSequence::length() const noexcept { return s.length(); }
inline size_t BasicNtSequence::max_size() const noexcept { return s.max_size(); }
inline void BasicNtSequence::resize (size_t n) { s.resize(n); }
inline void BasicNtSequence::resize (size_t n, char base) { Base::validate(base); s.resize(n, base); }
inline size_t BasicNtSequence::capacity() const noexcept { return s.capacity(); }
inline void BasicNtSequence::reserve(size_t n) { s.reserve(n); }
inline void BasicNtSequence::clear() noexcept { s.clear(); }
inline bool BasicNtSequence::empty() const noexcept { return s.empty(); }
inline void BasicNtSequence::shrink_to_fit() { s.shrink_to_fit(); }

inline BasicNtSequence::Base BasicNtSequence::operator[](size_t pos) { return Base(s[pos]); }
inline const char& BasicNtSequence::operator[](size_t pos) const { return s[pos]; }
inline BasicNtSequence::Base BasicNtSequence::at(size_t pos) { return Base(s.at(pos)); }
inline const char& BasicNtSequence::at(size_t pos) const { return s.at(pos); }
inline BasicNtSequence::Base BasicNtSequence::back() { return Base(s.back()); }
inline const char& BasicNtSequence::back() const { return s.back(); }
inline BasicNtSequence::Base BasicNtSequence::front() { return Base(s.front()); }
inline const char& BasicNtSequence::front() const { return s.front(); }

inline BasicNtSequence& BasicNtSequence::operator+=(const BasicNtSequence& rhs) { s += rhs.s; return *this; }
inline BasicNtSequence& BasicNtSequence::operator+=(const std::string& rhs) { validate(rhs); s += rhs; return *this; }
inline BasicNtSequence& BasicNtSequence::operator+=(const char* rhs) { validate(rhs); s += rhs; return *this; }
inline BasicNtSequence& BasicNtSequence::operator+=(char rhs) { Base::validate(rhs); s += rhs; return *this; }
inline BasicNtSequence& BasicNtSequence::operator+=(std::initializer_list<char> rhs) { validate(rhs); s += rhs; return *this; }

inline BasicNtSequence& BasicNtSequence::append(const BasicNtSequence& seq) { s.append(seq.s); return *this; }
inline BasicNtSequence& BasicNtSequence::append(const BasicNtSequence& seq, size_t subpos, size_t sublen) { s.append(seq.s, subpos, sublen); return *this; }
inline BasicNtSequence& BasicNtSequence::append(const std::string& bases) { validate(bases); s.append(bases); return *this; }
inline BasicNtSequence& BasicNtSequence::append(const std::string& bases, size_t subpos, size_t sublen) { validate(bases); s.append(bases, subpos, sublen); return *this; }
inline BasicNtSequence& BasicNtSequence::append(const char* bases) { validate(bases); s.append(bases); return *this; }
inline BasicNtSequence& BasicNtSequence::append(const char* bases, size_t n) { validate(bases); s.append(bases, n); return *this; }
inline BasicNtSequence& BasicNtSequence::append(size_t n, char base) { Base::validate(base); s.append(n, base); return *this; }
template <class InputIterator>
inline BasicNtSequence& BasicNtSequence::append(InputIterator first, InputIterator last) { s.append(first, last); validate(s); return *this; }
inline BasicNtSequence& BasicNtSequence::append(std::initializer_list<char> bases) { validate(bases); s.append(bases); return *this; }
inline void BasicNtSequence::push_back(char base) { Base::validate(base); s.push_back(base); }
inline BasicNtSequence& BasicNtSequence::assign(const BasicNtSequence& seq) { s.assign(seq.s); return *this; }
inline BasicNtSequence& BasicNtSequence::assign(const BasicNtSequence& seq, size_t subpos, size_t sublen) { s.assign(seq.s, subpos, sublen); return *this; }
inline BasicNtSequence& BasicNtSequence::assign(const std::string& bases) { validate(bases); s.assign(bases); return *this; }
inline BasicNtSequence& BasicNtSequence::assign(const std::string& bases, size_t subpos, size_t sublen) { validate(bases); s.assign(bases, subpos, npos); return *this; }
inline BasicNtSequence& BasicNtSequence::assign(const char* bases) { validate(bases); s.assign(bases); return *this; }
inline BasicNtSequence& BasicNtSequence::assign(const char* bases, size_t n) { validate(bases); s.assign(bases, n); return *this; }
inline BasicNtSequence& BasicNtSequence::assign(size_t n, char base) { Base::validate(base); s.assign(n, base); return *this; }
template <class InputIterator>
inline BasicNtSequence& BasicNtSequence::assign(InputIterator first, InputIterator last) { s.assign(first, last); validate(s); return *this; }
inline BasicNtSequence& BasicNtSequence::assign(std::initializer_list<char> bases) { validate(bases); s.assign(bases); return *this; }
inline BasicNtSequence& BasicNtSequence::assign(BasicNtSequence&& seq) noexcept { s.assign(seq.s); return *this; }
inline BasicNtSequence& BasicNtSequence::assign(std::string&& bases) noexcept { validate(bases); s.assign(bases); return *this; }

inline BasicNtSequence& BasicNtSequence::insert(size_t pos, const BasicNtSequence& seq) { s.insert(pos, seq.s); return *this; }
inline BasicNtSequence& BasicNtSequence::insert(size_t pos, const BasicNtSequence& seq, size_t subpos, size_t sublen) { s.insert(pos, seq.s, subpos, sublen); return *this; }
inline BasicNtSequence& BasicNtSequence::insert(size_t pos, const std::string& bases) { validate(bases); s.insert(pos, bases); return *this; }
inline BasicNtSequence& BasicNtSequence::insert(size_t pos, const std::string& bases, size_t subpos, size_t sublen) { validate(bases); s.insert(pos, bases, subpos, sublen); return *this; }
inline BasicNtSequence& BasicNtSequence::insert(size_t pos, const char* bases) { validate(bases); s.insert(pos, bases); return *this; }
inline BasicNtSequence& BasicNtSequence::insert(size_t pos, const char* bases, size_t n) { validate(bases); s.insert(pos,bases, n); return *this; }
inline BasicNtSequence& BasicNtSequence::insert(size_t pos, size_t n, char base) { Base::validate(base); s.insert(pos, n, base); return *this; }
inline BasicNtSequence::iterator BasicNtSequence::insert(const_iterator p, size_t n, char base) { Base::validate(base); return s.insert(p, n, base); }
inline BasicNtSequence::iterator BasicNtSequence::insert(const_iterator p, char base) { Base::validate(base); return s.insert(p, base); }
template <class InputIterator>
inline BasicNtSequence::iterator BasicNtSequence::insert(iterator p, InputIterator first, InputIterator last) { auto ret = s.insert(p, first, last); validate(s); return ret; }
inline BasicNtSequence& BasicNtSequence::insert(const_iterator p, std::initializer_list<char> bases) { validate(bases); s.insert(p, bases); return *this; }

inline BasicNtSequence& BasicNtSequence::erase (size_t pos, size_t len) { s.erase(pos, len); return *this; }
inline BasicNtSequence::iterator BasicNtSequence::erase (iterator p) { return s.erase(p); }
inline BasicNtSequence::iterator BasicNtSequence::erase (iterator first, iterator last) { return s.erase(first, last); }

inline BasicNtSequence& BasicNtSequence::replace(size_t pos, size_t len, const BasicNtSequence& seq) { s.replace(pos, len, seq.s); return *this; }
inline BasicNtSequence& BasicNtSequence::replace(size_t pos, size_t len, const std::string& bases) { validate(bases); s.replace(pos, len, bases); return *this; }
inline BasicNtSequence& BasicNtSequence::replace(const_iterator i1, const_iterator i2, const BasicNtSequence& seq) { s.replace(i1, i2, seq.s); return *this; }
inline BasicNtSequence& BasicNtSequence::replace(const_iterator i1, const_iterator i2, const std::string& bases) { validate(bases); s.replace(i1, i2, bases); return *this; }
inline BasicNtSequence& BasicNtSequence::replace(size_t pos, size_t len, const BasicNtSequence& seq, size_t subpos, size_t sublen) { s.replace(pos, len, seq.s, subpos, sublen); return *this; }
inline BasicNtSequence& BasicNtSequence::replace(size_t pos, size_t len, const std::string& bases, size_t subpos, size_t sublen) { validate(bases); s.replace(pos, len, bases, subpos, sublen); return *this; }
inline BasicNtSequence& BasicNtSequence::replace(size_t pos, size_t len, const char* bases) { validate(bases); s.replace(pos, len, bases); return *this; }
inline BasicNtSequence& BasicNtSequence::replace(const_iterator i1, const_iterator i2, const char* bases) { validate(bases); s.replace(i1, i2, bases); return *this; }
inline BasicNtSequence& BasicNtSequence::replace(size_t pos, size_t len, const char* bases, size_t n) { validate(bases); s.replace(pos, len, bases, n); return *this; }
inline BasicNtSequence& BasicNtSequence::replace(const_iterator i1, const_iterator i2, const char* bases, size_t n) { validate(bases); s.replace(i1, i2, bases, n); return *this; }
inline BasicNtSequence& BasicNtSequence::replace(size_t pos, size_t len, size_t n, char base) { Base::validate(base); s.replace(pos, len, n, base); return *this; }
inline BasicNtSequence& BasicNtSequence::replace(const_iterator i1, const_iterator i2, size_t n, char base) { Base::validate(base); s.replace(i1, i2, n, base); return *this; }
template <class InputIterator>
inline BasicNtSequence& BasicNtSequence::replace(const_iterator i1, const_iterator i2, InputIterator first, InputIterator last) { s.replace(i1, i2, first, last); validate(s); return *this; }
inline BasicNtSequence& BasicNtSequence::replace(const_iterator i1, const_iterator i2, std::initializer_list<char> bases) { validate(bases); s.replace(i1, i2, bases); return *this; }

inline void BasicNtSequence::swap(BasicNtSequence& seq) { s.swap(seq.s); }
inline void BasicNtSequence::swap(std::string& bases) { validate(bases); s.swap(bases); }

inline void BasicNtSequence::pop_back() { s.pop_back(); }

inline const char* BasicNtSequence::c_str() const noexcept { return s.c_str(); }
inline const char* BasicNtSequence::data() const noexcept { return s.data(); }

inline BasicNtSequence::allocator_type BasicNtSequence::get_allocator() const noexcept { return s.get_allocator(); }

inline size_t BasicNtSequence::copy(char* bases, size_t len, size_t pos) const { return s.copy(bases, len, pos); }

inline size_t BasicNtSequence::find(const BasicNtSequence& seq, size_t pos) const noexcept { return s.find(seq.s, pos); }
inline size_t BasicNtSequence::find(const std::string& bases, size_t pos) const noexcept { return s.find(bases, pos); }
inline size_t BasicNtSequence::find(const char* bases, size_t pos) const { return s.find(bases, pos); }
inline size_t BasicNtSequence::find(const char* bases, size_t pos, size_type n) const { return s.find(bases, pos, n); }
inline size_t BasicNtSequence::find(char base, size_t pos) const noexcept { return s.find(base, pos); }

inline size_t BasicNtSequence::rfind(const BasicNtSequence& seq, size_t pos) const noexcept { return s.rfind(seq.s, pos); }
inline size_t BasicNtSequence::rfind(const std::string& bases, size_t pos) const noexcept { return s.rfind(bases, pos); }
inline size_t BasicNtSequence::rfind(const char* bases, size_t pos) const { return s.rfind(bases, pos); }
inline size_t BasicNtSequence::rfind(const char* bases, size_t pos, size_t n) const { return s.rfind(bases, pos, n); }
inline size_t BasicNtSequence::rfind(char base, size_t pos) const noexcept { return s.rfind(base, pos); }

inline size_t BasicNtSequence::find_first_of(const BasicNtSequence& seq, size_t pos) const noexcept { return s.find_first_of(seq.s, pos); }
inline size_t BasicNtSequence::find_first_of(const std::string& bases, size_t pos) const noexcept { return s.find_first_of(bases, pos); }
inline size_t BasicNtSequence::find_first_of(const char* bases, size_t pos) const { return s.find_first_of(bases, pos); }
inline size_t BasicNtSequence::find_first_of(const char* bases, size_t pos, size_t n) const { return s.find_first_of(bases, pos, n); }
inline size_t BasicNtSequence::find_first_of(char base, size_t pos) const noexcept { return s.find_first_of(base, pos); }

inline size_t BasicNtSequence::find_last_of(const BasicNtSequence& seq, size_t pos) const noexcept { return s.find_last_of(seq.s, pos); }
inline size_t BasicNtSequence::find_last_of(const std::string& bases, size_t pos) const noexcept { return s.find_last_of(bases, pos); }
inline size_t BasicNtSequence::find_last_of(const char* bases, size_t pos) const { return s.find_last_of(bases, pos); }
inline size_t BasicNtSequence::find_last_of(const char* bases, size_t pos, size_t n) const { return s.find_last_of(bases, pos, n); }
inline size_t BasicNtSequence::find_last_of(char base, size_t pos) const noexcept { return s.find_last_of(base, pos); }

inline size_t BasicNtSequence::find_first_not_of(const BasicNtSequence& seq, size_t pos) const noexcept { return s.find_first_not_of(seq.s, pos); }
inline size_t BasicNtSequence::find_first_not_of(const std::string& bases, size_t pos) const noexcept { return s.find_first_not_of(bases, pos); }
inline size_t BasicNtSequence::find_first_not_of(const char* bases, size_t pos) const { return s.find_first_not_of(bases, pos); }
inline size_t BasicNtSequence::find_first_not_of(const char* bases, size_t pos, size_t n) const { return s.find_first_not_of(bases, pos, n); }
inline size_t BasicNtSequence::find_first_not_of(char base, size_t pos) const noexcept { return s.find_first_not_of(base, pos); }

inline size_t BasicNtSequence::find_last_not_of(const BasicNtSequence& seq, size_t pos) const noexcept { return s.find_last_not_of(seq.s, pos); }
inline size_t BasicNtSequence::find_last_not_of(const std::string& bases, size_t pos) const noexcept { return s.find_last_not_of(bases, pos); }
inline size_t BasicNtSequence::find_last_not_of(const char* bases, size_t pos) const { return s.find_last_not_of(bases, pos); }
inline size_t BasicNtSequence::find_last_not_of(const char* bases, size_t pos, size_t n) const { return s.find_last_not_of(bases, pos, n); }
inline size_t BasicNtSequence::find_last_not_of(char base, size_t pos) const noexcept { return s.find_last_not_of(base, pos); }

inline BasicNtSequence BasicNtSequence::substr(size_t pos, size_t len) const { BasicNtSequence seq; seq.set(s.substr(pos, len)); return seq; }

inline int BasicNtSequence::compare(const BasicNtSequence& seq) const noexcept { return s.compare(seq.s); }
inline int BasicNtSequence::compare(const std::string& bases) const noexcept { return s.compare(bases); }
inline int BasicNtSequence::compare(size_t pos, size_t len, const BasicNtSequence& seq) const { return s.compare(pos, len, seq.s); }
inline int BasicNtSequence::compare(size_t pos, size_t len, const std::string& bases) const { return s.compare(pos, len, bases); }
inline int BasicNtSequence::compare(size_t pos, size_t len, const BasicNtSequence& seq, size_t subpos, size_t sublen) const { return s.compare(pos, len, seq.s, subpos, sublen); }
inline int BasicNtSequence::compare(size_t pos, size_t len, const std::string& bases, size_t subpos, size_t sublen) const { return s.compare(pos, len, bases, subpos, sublen); }
inline int BasicNtSequence::compare(const char* bases) const { return s.compare(bases); }
inline int BasicNtSequence::compare(size_t pos, size_t len, const char* bases) const { return s.compare(pos, len, bases); }
inline int BasicNtSequence::compare(size_t pos, size_t len, const char* bases, size_t n) const { return s.compare(pos, len, bases, n); }

inline BasicNtSequence operator+(const BasicNtSequence& lhs, const BasicNtSequence& rhs) { return lhs.s + rhs.s; }
inline BasicNtSequence operator+(BasicNtSequence&& lhs, BasicNtSequence&& rhs) { return lhs.s + rhs.s; }
inline BasicNtSequence operator+(BasicNtSequence&& lhs, const BasicNtSequence& rhs) { return lhs.s + rhs.s; }
inline BasicNtSequence operator+(const BasicNtSequence& lhs, BasicNtSequence&& rhs) { return lhs.s + rhs.s; }
inline BasicNtSequence operator+(const std::string& lhs, const BasicNtSequence& rhs) { return lhs + rhs.s; }
inline BasicNtSequence operator+(const BasicNtSequence& lhs, const std::string& rhs) { return lhs.s + rhs; }
inline BasicNtSequence operator+(const std::string& lhs, BasicNtSequence&& rhs) { return lhs + rhs.s; }
inline BasicNtSequence operator+(const BasicNtSequence& lhs, std::string&& rhs) { return lhs.s + rhs; }
inline BasicNtSequence operator+(std::string&& lhs, const BasicNtSequence& rhs) { return lhs + rhs.s; }
inline BasicNtSequence operator+(BasicNtSequence&& lhs, const std::string& rhs) { return lhs.s + rhs; }
inline BasicNtSequence operator+(BasicNtSequence&& lhs, std::string&& rhs) { return lhs.s + rhs; }
inline BasicNtSequence operator+(std::string&& lhs, BasicNtSequence&& rhs) { return lhs + rhs.s; }
inline BasicNtSequence operator+(const char* lhs, const BasicNtSequence& rhs) { return lhs + rhs.s; }
inline BasicNtSequence operator+(const BasicNtSequence& lhs, const char* rhs) { return lhs.s + rhs; }
inline BasicNtSequence operator+(const char* lhs, BasicNtSequence&& rhs) { return lhs + rhs.s; }
inline BasicNtSequence operator+(BasicNtSequence&& lhs, const char* rhs) { return lhs.s + rhs; }
inline BasicNtSequence operator+(char lhs, const BasicNtSequence& rhs) { return lhs + rhs.s; }
inline BasicNtSequence operator+(const BasicNtSequence& lhs, char rhs) { return lhs.s + rhs; }
inline BasicNtSequence operator+(char lhs, BasicNtSequence&& rhs) { return lhs + rhs.s; }
inline BasicNtSequence operator+(BasicNtSequence&& lhs, char rhs) { return lhs.s + rhs; }

inline bool operator==(const BasicNtSequence& lhs, const BasicNtSequence& rhs) { return lhs.s == rhs.s; }
inline bool operator==(const std::string& lhs, const BasicNtSequence& rhs) { return lhs == rhs.s; }
inline bool operator==(const BasicNtSequence& lhs, const std::string& rhs) { return lhs.s == rhs; }
inline bool operator==(const char* lhs, const BasicNtSequence& rhs) { return lhs == rhs.s; }
inline bool operator==(const BasicNtSequence& lhs, const char* rhs) { return lhs.s == rhs; }

inline bool operator!=(const BasicNtSequence& lhs, const BasicNtSequence& rhs) { return lhs.s != rhs.s; }
inline bool operator!=(const std::string& lhs, const BasicNtSequence& rhs) { return lhs != rhs.s; }
inline bool operator!=(const BasicNtSequence& lhs, const std::string& rhs) { return lhs.s != rhs; }
inline bool operator!=(const char* lhs, const BasicNtSequence& rhs) { return lhs != rhs.s; }
inline bool operator!=(const BasicNtSequence& lhs, const char* rhs) { return lhs.s != rhs; }

inline bool operator<(const BasicNtSequence& lhs, const BasicNtSequence& rhs) { return lhs.s < rhs.s; }
inline bool operator<(const std::string& lhs, const BasicNtSequence& rhs) { return lhs < rhs.s; }
inline bool operator<(const BasicNtSequence& lhs, const std::string& rhs) { return lhs.s < rhs; }
inline bool operator<(const char* lhs, const BasicNtSequence& rhs) { return lhs < rhs.s; }
inline bool operator<(const BasicNtSequence& lhs, const char* rhs) { return lhs.s < rhs; }

inline bool operator<=(const BasicNtSequence& lhs, const BasicNtSequence& rhs) { return lhs.s <= rhs.s; }
inline bool operator<=(const std::string& lhs, const BasicNtSequence& rhs) { return lhs <= rhs.s; }
inline bool operator<=(const BasicNtSequence& lhs, const std::string& rhs) { return lhs.s <= rhs; }
inline bool operator<=(const char* lhs, const BasicNtSequence& rhs) { return lhs <= rhs.s; }
inline bool operator<=(const BasicNtSequence& lhs, const char* rhs) { return lhs.s <= rhs; }

inline bool operator>(const BasicNtSequence& lhs, const BasicNtSequence& rhs) { return lhs.s > rhs.s; }
inline bool operator>(const std::string& lhs, const BasicNtSequence& rhs) { return lhs > rhs.s; }
inline bool operator>(const BasicNtSequence& lhs, const std::string& rhs) { return lhs.s > rhs; }
inline bool operator>(const char* lhs, const BasicNtSequence& rhs) { return lhs > rhs.s; }
inline bool operator>(const BasicNtSequence& lhs, const char* rhs) { return lhs.s > rhs; }

inline bool operator>=(const BasicNtSequence& lhs, const BasicNtSequence& rhs) { return lhs.s >= rhs.s; }
inline bool operator>=(const std::string& lhs, const BasicNtSequence& rhs) { return lhs >= rhs.s; }
inline bool operator>=(const BasicNtSequence& lhs, const std::string& rhs) { return lhs.s >= rhs; }
inline bool operator>=(const char* lhs, const BasicNtSequence& rhs) { return lhs >= rhs.s; }
inline bool operator>=(const BasicNtSequence& lhs, const char* rhs) { return lhs.s >= rhs; }

inline void swap(BasicNtSequence& x, BasicNtSequence& y) { swap(x.s, y.s); }

inline std::istream& operator>>(std::istream& is, BasicNtSequence& rhs) { is >> rhs.s; return is; }
inline std::ostream& operator<<(std::ostream& os, const BasicNtSequence& rhs) { os << rhs.s; return os; }

inline std::istream& getline(std::istream& is, BasicNtSequence& seq, char delim) { auto& ret = getline(is, seq.s, delim); BasicNtSequence::validate(seq.s); return ret; }
inline std::istream& getline(std::istream&& is, BasicNtSequence& seq, char delim) { auto& ret = getline(is, seq.s, delim); BasicNtSequence::validate(seq.s); return ret; }
inline std::istream& getline(std::istream& is, BasicNtSequence& seq) { auto& ret = getline(is, seq.s); BasicNtSequence::validate(seq.s); return ret; }
inline std::istream& getline(std::istream&& is, BasicNtSequence& seq) { auto& ret = getline(is, seq.s); BasicNtSequence::validate(seq.s); return ret; }

inline void BasicNtSequence::validate(const std::string& bases) {
    std::for_each(bases.begin(), bases.end(),
        [] (char base) { Base::validate(base); }
    );
}

inline void BasicNtSequence::validate(const char* bases) {
    for (unsigned i = 0; bases[i] != '\0'; ++i) {
        Base::validate(bases[i]);  
    }
}

inline void BasicNtSequence::validate(std::initializer_list<char> bases) {
    std::for_each(bases.begin(), bases.end(),
        [] (char base) { Base::validate(base); }
    );
}

inline void BasicNtSequence::set(std::string bases) { s = std::move(bases); }
inline void BasicNtSequence::set(const char* bases) { s = std::string(bases); }

} // namespace btl

#endif
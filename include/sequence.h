#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>

namespace btl {

class Sequence {

public:

    struct Base {

        Base(char& base);
        Base(const Base& base);

        Base& operator=(char base);
        Base& operator=(const Base& base);

        operator char() const;

        void complement();
        char getComplement() const;
        char operator~() const;

        void validate();
        void capitalize();

        static inline void validate(char base);
        static inline char capitalize(char base);

        friend class Sequence;

        char& b;
        static const char COMPLEMENTS[256];
        static const char CAPITALS[256];
    };

    static inline const size_t npos = std::string::npos;

    Sequence();
    Sequence(const Sequence& seq);
    Sequence(std::string bases);
    Sequence(const Sequence& seq, size_t pos, size_t len = npos);
    Sequence(const std::string& bases, size_t pos, size_t len = npos);
    Sequence(const char* bases);
    Sequence(const char* bases, size_t n);
    Sequence(size_t n, char base);
    template <class InputIterator>
    Sequence(InputIterator first, InputIterator last);
    Sequence(std::initializer_list<char> bases);
    Sequence(Sequence&& seq) noexcept;

    virtual ~Sequence() {}

    Sequence& operator=(const Sequence& seq);
    Sequence& operator=(std::string bases);
    Sequence& operator=(const char* bases);
    Sequence& operator=(std::initializer_list<char> bases);
    Sequence& operator=(Sequence&& seq) noexcept;
    Sequence& operator=(std::string&& bases);

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
    char operator[](size_t pos) const;
    Base at(size_t pos);
    char at(size_t pos) const;
    Base back();
    char back() const;
    Base front();
    char front() const;

    Sequence& operator+=(const Sequence& rhs);
    Sequence& operator+=(const std::string& rhs);
    Sequence& operator+=(const char* rhs);
    Sequence& operator+=(char rhs);
    Sequence& operator+=(std::initializer_list<char> rhs);

    Sequence& append(const Sequence& seq);
    Sequence& append(const Sequence& seq, size_t subpos, size_t sublen);
    Sequence& append(const std::string& bases);
    Sequence& append(const std::string& bases, size_t subpos, size_t sublen);
    Sequence& append(const char* bases);
    Sequence& append(const char* bases, size_t n);
    Sequence& append(size_t n, char base);
    template <class InputIterator>
    Sequence& append(InputIterator first, InputIterator last);
    Sequence& append(std::initializer_list<char> bases);

    void push_back(char base);

    Sequence& assign(const Sequence& seq);
    Sequence& assign(const Sequence& seq, size_t subpos, size_t sublen = npos);
    Sequence& assign(const std::string& bases);
    Sequence& assign(const std::string& bases, size_t subpos, size_t sublen = npos);
    Sequence& assign(const char* bases);
    Sequence& assign(const char* bases, size_t n);
    Sequence& assign(size_t n, char base);
    template <class InputIterator>
    Sequence& assign(InputIterator first, InputIterator last);
    Sequence& assign(std::initializer_list<char> bases);
    Sequence& assign(Sequence&& seq) noexcept;
    Sequence& assign(std::string&& bases) noexcept;

    Sequence& insert(size_t pos, const Sequence& seq);
    Sequence& insert(size_t pos, const Sequence& seq, size_t subpos, size_t sublen = npos);
    Sequence& insert(size_t pos, const std::string& bases);
    Sequence& insert(size_t pos, const std::string& bases, size_t subpos, size_t sublen = npos);
    Sequence& insert(size_t pos, const char* bases);
    Sequence& insert(size_t pos, const char* bases, size_t n);
    Sequence& insert(size_t pos, size_t n, char base);
    iterator insert(const_iterator p, size_t n, char base);
    iterator insert(const_iterator p, char base);
    template <class InputIterator>
    iterator insert(iterator p, InputIterator first, InputIterator last);
    Sequence& insert(const_iterator p, std::initializer_list<char> bases);

    Sequence& erase(size_t pos = 0, size_t len = npos);
    iterator erase(iterator p);
    iterator erase(iterator first, iterator last);

    Sequence& replace(size_t pos, size_t len, const Sequence& seq);
    Sequence& replace(size_t pos, size_t len, const std::string& bases);
    Sequence& replace(const_iterator i1, const_iterator i2, const Sequence& seq);
    Sequence& replace(const_iterator i1, const_iterator i2, const std::string& bases);
    Sequence& replace(size_t pos, size_t len, const Sequence& seq, size_t subpos, size_t sublen = npos);
    Sequence& replace(size_t pos, size_t len, const std::string& bases, size_t subpos, size_t sublen = npos);
    Sequence& replace(size_t pos, size_t len, const char* bases);
    Sequence& replace(const_iterator i1, const_iterator i2, const char* bases);
    Sequence& replace(size_t pos, size_t len, const char* bases, size_t n);
    Sequence& replace(const_iterator i1, const_iterator i2, const char* bases, size_t n);
    Sequence& replace(size_t pos, size_t len, size_t n, char base);
    Sequence& replace(const_iterator i1, const_iterator i2, size_t n, char base);
    template <class InputIterator>
    Sequence& replace(const_iterator i1, const_iterator i2, InputIterator first, InputIterator last);
    Sequence& replace(const_iterator i1, const_iterator i2, std::initializer_list<char> bases);

    void swap(Sequence& seq);
    void swap(std::string& bases);

    void pop_back();

    const char* c_str() const noexcept;
    const char* data() const noexcept;

    using allocator_type = std::string::allocator_type;

    allocator_type get_allocator() const noexcept;

    size_t copy(char* bases, size_t len, size_t pos = 0) const;

    using size_type = std::string::size_type;

    size_t find(const Sequence& seq, size_t pos = 0) const noexcept;
    size_t find(const std::string& bases, size_t pos = 0) const noexcept;
    size_t find(const char* bases, size_t pos = 0) const;
    size_t find(const char* bases, size_t pos, size_type n) const;
    size_t find(char base, size_t pos = 0) const noexcept;

    size_t rfind(const Sequence& seq, size_t pos = npos) const noexcept;
    size_t rfind(const std::string& bases, size_t pos = npos) const noexcept;
    size_t rfind(const char* bases, size_t pos = npos) const;
    size_t rfind(const char* bases, size_t pos, size_t n) const;
    size_t rfind(char base, size_t pos = npos) const noexcept;

    size_t find_first_of(const Sequence& seq, size_t pos = 0) const noexcept;
    size_t find_first_of(const std::string& bases, size_t pos = 0) const noexcept;
    size_t find_first_of(const char* bases, size_t pos = 0) const;
    size_t find_first_of(const char* bases, size_t pos, size_t n) const;
    size_t find_first_of(char base, size_t pos = 0) const noexcept;

    size_t find_last_of(const Sequence& seq, size_t pos = npos) const noexcept;
    size_t find_last_of(const std::string& bases, size_t pos = npos) const noexcept;
    size_t find_last_of(const char* bases, size_t pos = npos) const;
    size_t find_last_of(const char* bases, size_t pos, size_t n) const;
    size_t find_last_of(char base, size_t pos = npos) const noexcept;

    size_t find_first_not_of(const Sequence& seq, size_t pos = 0) const noexcept;
    size_t find_first_not_of(const std::string& bases, size_t pos = 0) const noexcept;
    size_t find_first_not_of(const char* bases, size_t pos = 0) const;
    size_t find_first_not_of(const char* bases, size_t pos, size_t n) const;
    size_t find_first_not_of(char base, size_t pos = 0) const noexcept;

    size_t find_last_not_of(const Sequence& seq, size_t pos = npos) const noexcept;
    size_t find_last_not_of(const std::string& bases, size_t pos = npos) const noexcept;
    size_t find_last_not_of(const char* bases, size_t pos = npos) const;
    size_t find_last_not_of(const char* bases, size_t pos, size_t n) const;
    size_t find_last_not_of(char base, size_t pos = npos) const noexcept;

    Sequence substr(size_t pos = 0, size_t len = npos) const;

    int compare(const Sequence& seq) const noexcept;
    int compare(const std::string& bases) const noexcept;
    int compare(size_t pos, size_t len, const Sequence& seq) const;
    int compare(size_t pos, size_t len, const std::string& bases) const;
    int compare(size_t pos, size_t len, const Sequence& seq, size_t subpos, size_t sublen = npos) const;
    int compare(size_t pos, size_t len, const std::string& bases, size_t subpos, size_t sublen = npos) const;
    int compare(const char* bases) const;
    int compare(size_t pos, size_t len, const char* bases) const;
    int compare(size_t pos, size_t len, const char* bases, size_t n) const;

    friend Sequence operator+(const Sequence& lhs, const Sequence& rhs);
    friend Sequence operator+(Sequence&& lhs, Sequence&& rhs);
    friend Sequence operator+(Sequence&& lhs, const Sequence& rhs);
    friend Sequence operator+(const Sequence& lhs, Sequence&& rhs);
    friend Sequence operator+(const std::string& lhs, const Sequence& rhs);
    friend Sequence operator+(const Sequence& lhs, const std::string& rhs);
    friend Sequence operator+(const std::string& lhs, Sequence&& rhs);
    friend Sequence operator+(const Sequence& lhs, std::string&& rhs);
    friend Sequence operator+(std::string&& lhs, const Sequence& rhs);
    friend Sequence operator+(Sequence&& lhs, const std::string& rhs);
    friend Sequence operator+(Sequence&& lhs, std::string&& rhs);
    friend Sequence operator+(std::string&& lhs, Sequence&& rhs);
    friend Sequence operator+(const char* lhs, const Sequence& rhs);
    friend Sequence operator+(const Sequence& lhs, const char* rhs);
    friend Sequence operator+(const char* lhs, Sequence&& rhs);
    friend Sequence operator+(Sequence&& lhs, const char* rhs);
    friend Sequence operator+(char lhs, const Sequence& rhs);
    friend Sequence operator+(const Sequence& lhs, char rhs);
    friend Sequence operator+(char lhs, Sequence&& rhs);
    friend Sequence operator+(Sequence&& lhs, char rhs);

    friend bool operator==(const Sequence& lhs, const Sequence& rhs);
    friend bool operator==(const std::string& lhs, const Sequence& rhs);
    friend bool operator==(const Sequence& lhs, const std::string& rhs);
    friend bool operator==(const char* lhs, const Sequence& rhs);
    friend bool operator==(const Sequence& lhs, const char* rhs);

    friend bool operator!=(const Sequence& lhs, const Sequence& rhs);
    friend bool operator!=(const std::string& lhs, const Sequence& rhs);
    friend bool operator!=(const Sequence& lhs, const std::string& rhs);
    friend bool operator!=(const char* lhs, const Sequence& rhs);
    friend bool operator!=(const Sequence& lhs, const char* rhs);

    friend bool operator<(const Sequence& lhs, const Sequence& rhs);
    friend bool operator<(const std::string& lhs, const Sequence& rhs);
    friend bool operator<(const Sequence& lhs, const std::string& rhs);
    friend bool operator<(const char* lhs, const Sequence& rhs);
    friend bool operator<(const Sequence& lhs, const char* rhs);

    friend bool operator<=(const Sequence& lhs, const Sequence& rhs);
    friend bool operator<=(const std::string& lhs, const Sequence& rhs);
    friend bool operator<=(const Sequence& lhs, const std::string& rhs);
    friend bool operator<=(const char* lhs, const Sequence& rhs);
    friend bool operator<=(const Sequence& lhs, const char* rhs);

    friend bool operator>(const Sequence& lhs, const Sequence& rhs);
    friend bool operator>(const std::string& lhs, const Sequence& rhs);
    friend bool operator>(const Sequence& lhs, const std::string& rhs);
    friend bool operator>(const char* lhs, const Sequence& rhs);
    friend bool operator>(const Sequence& lhs, const char* rhs);

    friend bool operator>=(const Sequence& lhs, const Sequence& rhs);
    friend bool operator>=(const std::string& lhs, const Sequence& rhs);
    friend bool operator>=(const Sequence& lhs, const std::string& rhs);
    friend bool operator>=(const char* lhs, const Sequence& rhs);
    friend bool operator>=(const Sequence& lhs, const char* rhs);

    friend void swap(Sequence& x, Sequence& y);

    friend std::istream& operator>>(std::istream& is, Sequence& rhs);
    friend std::ostream& operator<<(std::ostream& os, const Sequence& rhs);

    friend std::istream& getline(std::istream& is, Sequence& seq, char delim);
    friend std::istream& getline(std::istream&& is, Sequence& seq, char delim);
    friend std::istream& getline(std::istream& is, Sequence& seq);
    friend std::istream& getline(std::istream&& is, Sequence& seq);

    void reverseComplement();
    Sequence getReverseComplement();
    Sequence operator~();

private:

    static void validate(const std::string& bases);
    static void validate(const char* bases);
    static void validate(std::initializer_list<char> bases);

    void validate();
    void capitalize();

    // No validation
    void set(std::string bases);
    void set(const char* bases);

    std::string s;

};

inline Sequence::Base::Base(char& base): b(base) {}
inline Sequence::Base::Base(const Base& base) = default;

inline Sequence::Base& Sequence::Base::operator=(char base) { validate(base); b = base; return *this; }
inline Sequence::Base& Sequence::Base::operator=(const Base& base) { b = base.b; return *this; }

inline Sequence::Base::operator char() const { return b; }

inline void Sequence::Base::complement() { b = COMPLEMENTS[(unsigned char)b]; }
inline char Sequence::Base::getComplement() const { return COMPLEMENTS[(unsigned char)b]; }
inline char Sequence::Base::operator~() const { return getComplement(); }

inline void Sequence::Base::validate() { validate(b); }
inline void Sequence::Base::capitalize() { b = CAPITALS[(unsigned char)b]; }

inline void Sequence::Base::validate(char base) { assert(COMPLEMENTS[(unsigned char)base]); }
inline char Sequence::Base::capitalize(char base) { return CAPITALS[(unsigned char)base]; }

const inline char Sequence::Base::COMPLEMENTS[256] = {
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//       !    "    #    $    %    &    '    (    )    *    +    ,    -    .    /
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , '-' , '.',  0 ,

//  0    1    2    3    4    5    6    7    8    9    :    ;    <    =    >    ?
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//  @    A    B    C    D    E    F    G    H    I    J    K    L    M    N    O
    0 , 'T', 'V', 'G', 'H',  0 ,  0 , 'C', 'D',  0 ,  0 , 'M',  0 , 'K', 'N',  0 , 

//  P    Q    R    S    T    U    V    W    X    Y    Z    [    \    ]    ^    _
    0 ,  0 , 'Y', 'S', 'A', 'U', 'B', 'W',  0 , 'R',  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//  `    a    b    c    d    e    f    g    h    i    j    k    l    m    n    o
    0 , 't', 'v', 'g', 'h',  0 ,  0 , 'c', 'd',  0 ,  0 , 'm',  0 , 'k', 'n',  0 ,

//  p    q    r    s    t    u    v    w    x    y    z    {    |    }    ~   DEL
    0 ,  0 , 'y', 's', 'a', 'u', 'b', 'w',  0 , 'r',  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0
};

const inline char Sequence::Base::CAPITALS[256] = {
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//       !    "    #    $    %    &    '    (    )    *    +    ,    -    .    /
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , '-' , '.',  0 ,

//  0    1    2    3    4    5    6    7    8    9    :    ;    <    =    >    ?
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//  @    A    B    C    D    E    F    G    H    I    J    K    L    M    N    O
    0 , 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 

//  P    Q    R    S    T    U    V    W    X    Y    Z    [    \    ]    ^    _
   'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',  0 ,  0 ,  0 ,  0 ,  0 ,

//  `    a    b    c    d    e    f    g    h    i    j    k    l    m    n    o
    0 , 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',

//  p    q    r    s    t    u    v    w    x    y    z    {    |    }    ~   DEL
   'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',  0 ,  0 ,  0 ,  0 ,  0 ,

    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0
};

inline Sequence::Sequence() = default;
inline Sequence::Sequence(const Sequence& seq) = default;
inline Sequence::Sequence(std::string bases): s(std::move(bases)) { validate(); capitalize(); }
inline Sequence::Sequence(const Sequence& seq, size_t pos, size_t len): s(seq.s, pos, len) {}
inline Sequence::Sequence(const std::string& bases, size_t pos, size_t len): s(bases, pos, len) { validate(); capitalize(); }
inline Sequence::Sequence(const char* bases): s(bases) { validate(); capitalize(); }
inline Sequence::Sequence(const char* bases, size_t n): s(bases, n) { validate(); capitalize(); }
inline Sequence::Sequence(size_t n, char base): s(n, base) { validate(); capitalize(); }
template <class InputIterator>
inline Sequence::Sequence(InputIterator first, InputIterator last): s(first, last) { validate(); capitalize(); }
inline Sequence::Sequence(std::initializer_list<char> bases): s(bases) { validate(); capitalize(); }
inline Sequence::Sequence(Sequence&& seq) noexcept: s(std::move(seq.s)) {}

inline Sequence& Sequence::operator=(const Sequence& seq) = default;
inline Sequence& Sequence::operator=(std::string bases) { validate(bases); s = std::move(bases); capitalize(); return *this; }
inline Sequence& Sequence::operator=(const char* bases) { validate(bases); s = bases; capitalize(); return *this; }
inline Sequence& Sequence::operator=(std::initializer_list<char> bases) { validate(bases); s = bases; capitalize(); return *this; }
inline Sequence& Sequence::operator=(Sequence&& seq) noexcept { s = std::move(seq.s); return *this; }
inline Sequence& Sequence::operator=(std::string&& bases) { validate(bases); s = bases; capitalize(); return *this; }

inline Sequence::operator const std::string&() const noexcept { return s; }
inline Sequence::operator const char*() const noexcept { return s.c_str(); }

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
inline void Sequence::resize (size_t n, char base) { Base::validate(base); s.resize(n, Base::capitalize(base)); }
inline size_t Sequence::capacity() const noexcept { return s.capacity(); }
inline void Sequence::reserve(size_t n) { s.reserve(n); }
inline void Sequence::clear() noexcept { s.clear(); }
inline bool Sequence::empty() const noexcept { return s.empty(); }
inline void Sequence::shrink_to_fit() { s.shrink_to_fit(); }

inline Sequence::Base Sequence::operator[](size_t pos) { return Base(s[pos]); }
inline char Sequence::operator[](size_t pos) const { return s[pos]; }
inline Sequence::Base Sequence::at(size_t pos) { return Base(s.at(pos)); }
inline char Sequence::at(size_t pos) const { return s.at(pos); }
inline Sequence::Base Sequence::back() { return Base(s.back()); }
inline char Sequence::back() const { return s.back(); }
inline Sequence::Base Sequence::front() { return Base(s.front()); }
inline char Sequence::front() const { return s.front(); }

inline Sequence& Sequence::operator+=(const Sequence& rhs) { s += rhs.s; return *this; }
inline Sequence& Sequence::operator+=(const std::string& rhs) { validate(rhs); s += rhs; capitalize(); return *this; }
inline Sequence& Sequence::operator+=(const char* rhs) { validate(rhs); s += rhs; capitalize(); return *this; }
inline Sequence& Sequence::operator+=(char rhs) { Base::validate(rhs); s += Base::capitalize(rhs); return *this; }
inline Sequence& Sequence::operator+=(std::initializer_list<char> rhs) { validate(rhs); s += rhs; capitalize(); return *this; }

inline Sequence& Sequence::append(const Sequence& seq) { s.append(seq.s); return *this; }
inline Sequence& Sequence::append(const Sequence& seq, size_t subpos, size_t sublen) { s.append(seq.s, subpos, sublen); return *this; }
inline Sequence& Sequence::append(const std::string& bases) { validate(bases); s.append(bases); capitalize(); return *this; }
inline Sequence& Sequence::append(const std::string& bases, size_t subpos, size_t sublen) { validate(bases); s.append(bases, subpos, sublen); capitalize(); return *this; }
inline Sequence& Sequence::append(const char* bases) { validate(bases); s.append(bases); capitalize(); return *this; }
inline Sequence& Sequence::append(const char* bases, size_t n) { validate(bases); s.append(bases, n); capitalize(); return *this; }
inline Sequence& Sequence::append(size_t n, char base) { Base::validate(base); s.append(n, Base::capitalize(base)); return *this; }
template <class InputIterator>
inline Sequence& Sequence::append(InputIterator first, InputIterator last) { s.append(first, last); validate(); capitalize(); return *this; }
inline Sequence& Sequence::append(std::initializer_list<char> bases) { validate(bases); s.append(bases); capitalize(); return *this; }
inline void Sequence::push_back(char base) { Base::validate(base); s.push_back(Base::capitalize(base)); }
inline Sequence& Sequence::assign(const Sequence& seq) { s.assign(seq.s); return *this; }
inline Sequence& Sequence::assign(const Sequence& seq, size_t subpos, size_t sublen) { s.assign(seq.s, subpos, sublen); return *this; }
inline Sequence& Sequence::assign(const std::string& bases) { validate(bases); s.assign(bases); capitalize(); return *this; }
inline Sequence& Sequence::assign(const std::string& bases, size_t subpos, size_t sublen) { validate(bases); s.assign(bases, subpos, sublen); capitalize(); return *this; }
inline Sequence& Sequence::assign(const char* bases) { validate(bases); s.assign(bases); capitalize(); return *this; }
inline Sequence& Sequence::assign(const char* bases, size_t n) { validate(bases); s.assign(bases, n); capitalize(); return *this; }
inline Sequence& Sequence::assign(size_t n, char base) { Base::validate(base); s.assign(n, Base::capitalize(base)); return *this; }
template <class InputIterator>
inline Sequence& Sequence::assign(InputIterator first, InputIterator last) { s.assign(first, last); validate(); capitalize(); return *this; }
inline Sequence& Sequence::assign(std::initializer_list<char> bases) { validate(bases); s.assign(bases); capitalize(); return *this; }
inline Sequence& Sequence::assign(Sequence&& seq) noexcept { s.assign(seq.s); return *this; }
inline Sequence& Sequence::assign(std::string&& bases) noexcept { validate(bases); s.assign(bases); capitalize(); return *this; }

inline Sequence& Sequence::insert(size_t pos, const Sequence& seq) { s.insert(pos, seq.s); return *this; }
inline Sequence& Sequence::insert(size_t pos, const Sequence& seq, size_t subpos, size_t sublen) { s.insert(pos, seq.s, subpos, sublen); return *this; }
inline Sequence& Sequence::insert(size_t pos, const std::string& bases) { validate(bases); s.insert(pos, bases); capitalize(); return *this; }
inline Sequence& Sequence::insert(size_t pos, const std::string& bases, size_t subpos, size_t sublen) { validate(bases); s.insert(pos, bases, subpos, sublen); capitalize(); return *this; }
inline Sequence& Sequence::insert(size_t pos, const char* bases) { validate(bases); s.insert(pos, bases); capitalize(); return *this; }
inline Sequence& Sequence::insert(size_t pos, const char* bases, size_t n) { validate(bases); s.insert(pos,bases, n); capitalize(); return *this; }
inline Sequence& Sequence::insert(size_t pos, size_t n, char base) { Base::validate(base); s.insert(pos, n, Base::capitalize(base)); return *this; }
inline Sequence::iterator Sequence::insert(const_iterator p, size_t n, char base) { Base::validate(base); return s.insert(p, n, Base::capitalize(base)); }
inline Sequence::iterator Sequence::insert(const_iterator p, char base) { Base::validate(base); return s.insert(p, Base::capitalize(base)); }
template <class InputIterator>
inline Sequence::iterator Sequence::insert(iterator p, InputIterator first, InputIterator last) { auto ret = s.insert(p, first, last); validate(); capitalize(); return ret; }
inline Sequence& Sequence::insert(const_iterator p, std::initializer_list<char> bases) { validate(bases); s.insert(p, bases); capitalize(); return *this; }

inline Sequence& Sequence::erase (size_t pos, size_t len) { s.erase(pos, len); return *this; }
inline Sequence::iterator Sequence::erase (iterator p) { return s.erase(p); }
inline Sequence::iterator Sequence::erase (iterator first, iterator last) { return s.erase(first, last); }

inline Sequence& Sequence::replace(size_t pos, size_t len, const Sequence& seq) { s.replace(pos, len, seq.s); return *this; }
inline Sequence& Sequence::replace(size_t pos, size_t len, const std::string& bases) { validate(bases); s.replace(pos, len, bases); capitalize(); return *this; }
inline Sequence& Sequence::replace(const_iterator i1, const_iterator i2, const Sequence& seq) { s.replace(i1, i2, seq.s); return *this; }
inline Sequence& Sequence::replace(const_iterator i1, const_iterator i2, const std::string& bases) { validate(bases); s.replace(i1, i2, bases); capitalize(); return *this; }
inline Sequence& Sequence::replace(size_t pos, size_t len, const Sequence& seq, size_t subpos, size_t sublen) { s.replace(pos, len, seq.s, subpos, sublen); return *this; }
inline Sequence& Sequence::replace(size_t pos, size_t len, const std::string& bases, size_t subpos, size_t sublen) { validate(bases); s.replace(pos, len, bases, subpos, sublen); capitalize(); return *this; }
inline Sequence& Sequence::replace(size_t pos, size_t len, const char* bases) { validate(bases); s.replace(pos, len, bases); capitalize(); return *this; }
inline Sequence& Sequence::replace(const_iterator i1, const_iterator i2, const char* bases) { validate(bases); s.replace(i1, i2, bases); capitalize(); return *this; }
inline Sequence& Sequence::replace(size_t pos, size_t len, const char* bases, size_t n) { validate(bases); s.replace(pos, len, bases, n); capitalize(); return *this; }
inline Sequence& Sequence::replace(const_iterator i1, const_iterator i2, const char* bases, size_t n) { validate(bases); s.replace(i1, i2, bases, n); capitalize(); return *this; }
inline Sequence& Sequence::replace(size_t pos, size_t len, size_t n, char base) { Base::validate(base); s.replace(pos, len, n, Base::capitalize(base)); return *this; }
inline Sequence& Sequence::replace(const_iterator i1, const_iterator i2, size_t n, char base) { Base::validate(base); s.replace(i1, i2, n, Base::capitalize(base)); return *this; }
template <class InputIterator>
inline Sequence& Sequence::replace(const_iterator i1, const_iterator i2, InputIterator first, InputIterator last) { s.replace(i1, i2, first, last); validate(); capitalize(); return *this; }
inline Sequence& Sequence::replace(const_iterator i1, const_iterator i2, std::initializer_list<char> bases) { validate(bases); s.replace(i1, i2, bases); capitalize(); return *this; }

inline void Sequence::swap(Sequence& seq) { s.swap(seq.s); }
inline void Sequence::swap(std::string& bases) { validate(bases); s.swap(bases); capitalize(); }

inline void Sequence::pop_back() { s.pop_back(); }

inline const char* Sequence::c_str() const noexcept { return s.c_str(); }
inline const char* Sequence::data() const noexcept { return s.data(); }

inline Sequence::allocator_type Sequence::get_allocator() const noexcept { return s.get_allocator(); }

inline size_t Sequence::copy(char* bases, size_t len, size_t pos) const { return s.copy(bases, len, pos); }

inline size_t Sequence::find(const Sequence& seq, size_t pos) const noexcept { return s.find(seq.s, pos); }
inline size_t Sequence::find(const std::string& bases, size_t pos) const noexcept { return s.find(bases, pos); }
inline size_t Sequence::find(const char* bases, size_t pos) const { return s.find(bases, pos); }
inline size_t Sequence::find(const char* bases, size_t pos, size_type n) const { return s.find(bases, pos, n); }
inline size_t Sequence::find(char base, size_t pos) const noexcept { return s.find(base, pos); }

inline size_t Sequence::rfind(const Sequence& seq, size_t pos) const noexcept { return s.rfind(seq.s, pos); }
inline size_t Sequence::rfind(const std::string& bases, size_t pos) const noexcept { return s.rfind(bases, pos); }
inline size_t Sequence::rfind(const char* bases, size_t pos) const { return s.rfind(bases, pos); }
inline size_t Sequence::rfind(const char* bases, size_t pos, size_t n) const { return s.rfind(bases, pos, n); }
inline size_t Sequence::rfind(char base, size_t pos) const noexcept { return s.rfind(base, pos); }

inline size_t Sequence::find_first_of(const Sequence& seq, size_t pos) const noexcept { return s.find_first_of(seq.s, pos); }
inline size_t Sequence::find_first_of(const std::string& bases, size_t pos) const noexcept { return s.find_first_of(bases, pos); }
inline size_t Sequence::find_first_of(const char* bases, size_t pos) const { return s.find_first_of(bases, pos); }
inline size_t Sequence::find_first_of(const char* bases, size_t pos, size_t n) const { return s.find_first_of(bases, pos, n); }
inline size_t Sequence::find_first_of(char base, size_t pos) const noexcept { return s.find_first_of(base, pos); }

inline size_t Sequence::find_last_of(const Sequence& seq, size_t pos) const noexcept { return s.find_last_of(seq.s, pos); }
inline size_t Sequence::find_last_of(const std::string& bases, size_t pos) const noexcept { return s.find_last_of(bases, pos); }
inline size_t Sequence::find_last_of(const char* bases, size_t pos) const { return s.find_last_of(bases, pos); }
inline size_t Sequence::find_last_of(const char* bases, size_t pos, size_t n) const { return s.find_last_of(bases, pos, n); }
inline size_t Sequence::find_last_of(char base, size_t pos) const noexcept { return s.find_last_of(base, pos); }

inline size_t Sequence::find_first_not_of(const Sequence& seq, size_t pos) const noexcept { return s.find_first_not_of(seq.s, pos); }
inline size_t Sequence::find_first_not_of(const std::string& bases, size_t pos) const noexcept { return s.find_first_not_of(bases, pos); }
inline size_t Sequence::find_first_not_of(const char* bases, size_t pos) const { return s.find_first_not_of(bases, pos); }
inline size_t Sequence::find_first_not_of(const char* bases, size_t pos, size_t n) const { return s.find_first_not_of(bases, pos, n); }
inline size_t Sequence::find_first_not_of(char base, size_t pos) const noexcept { return s.find_first_not_of(base, pos); }

inline size_t Sequence::find_last_not_of(const Sequence& seq, size_t pos) const noexcept { return s.find_last_not_of(seq.s, pos); }
inline size_t Sequence::find_last_not_of(const std::string& bases, size_t pos) const noexcept { return s.find_last_not_of(bases, pos); }
inline size_t Sequence::find_last_not_of(const char* bases, size_t pos) const { return s.find_last_not_of(bases, pos); }
inline size_t Sequence::find_last_not_of(const char* bases, size_t pos, size_t n) const { return s.find_last_not_of(bases, pos, n); }
inline size_t Sequence::find_last_not_of(char base, size_t pos) const noexcept { return s.find_last_not_of(base, pos); }

inline Sequence Sequence::substr(size_t pos, size_t len) const { Sequence seq; seq.set(s.substr(pos, len)); return seq; }

inline int Sequence::compare(const Sequence& seq) const noexcept { return s.compare(seq.s); }
inline int Sequence::compare(const std::string& bases) const noexcept { return s.compare(bases); }
inline int Sequence::compare(size_t pos, size_t len, const Sequence& seq) const { return s.compare(pos, len, seq.s); }
inline int Sequence::compare(size_t pos, size_t len, const std::string& bases) const { return s.compare(pos, len, bases); }
inline int Sequence::compare(size_t pos, size_t len, const Sequence& seq, size_t subpos, size_t sublen) const { return s.compare(pos, len, seq.s, subpos, sublen); }
inline int Sequence::compare(size_t pos, size_t len, const std::string& bases, size_t subpos, size_t sublen) const { return s.compare(pos, len, bases, subpos, sublen); }
inline int Sequence::compare(const char* bases) const { return s.compare(bases); }
inline int Sequence::compare(size_t pos, size_t len, const char* bases) const { return s.compare(pos, len, bases); }
inline int Sequence::compare(size_t pos, size_t len, const char* bases, size_t n) const { return s.compare(pos, len, bases, n); }

inline Sequence operator+(const Sequence& lhs, const Sequence& rhs) { return lhs.s + rhs.s; }
inline Sequence operator+(Sequence&& lhs, Sequence&& rhs) { return lhs.s + rhs.s; }
inline Sequence operator+(Sequence&& lhs, const Sequence& rhs) { return lhs.s + rhs.s; }
inline Sequence operator+(const Sequence& lhs, Sequence&& rhs) { return lhs.s + rhs.s; }
inline Sequence operator+(const std::string& lhs, const Sequence& rhs) { return lhs + rhs.s; }
inline Sequence operator+(const Sequence& lhs, const std::string& rhs) { return lhs.s + rhs; }
inline Sequence operator+(const std::string& lhs, Sequence&& rhs) { return lhs + rhs.s; }
inline Sequence operator+(const Sequence& lhs, std::string&& rhs) { return lhs.s + rhs; }
inline Sequence operator+(std::string&& lhs, const Sequence& rhs) { return lhs + rhs.s; }
inline Sequence operator+(Sequence&& lhs, const std::string& rhs) { return lhs.s + rhs; }
inline Sequence operator+(Sequence&& lhs, std::string&& rhs) { return lhs.s + rhs; }
inline Sequence operator+(std::string&& lhs, Sequence&& rhs) { return lhs + rhs.s; }
inline Sequence operator+(const char* lhs, const Sequence& rhs) { return lhs + rhs.s; }
inline Sequence operator+(const Sequence& lhs, const char* rhs) { return lhs.s + rhs; }
inline Sequence operator+(const char* lhs, Sequence&& rhs) { return lhs + rhs.s; }
inline Sequence operator+(Sequence&& lhs, const char* rhs) { return lhs.s + rhs; }
inline Sequence operator+(char lhs, const Sequence& rhs) { return lhs + rhs.s; }
inline Sequence operator+(const Sequence& lhs, char rhs) { return lhs.s + rhs; }
inline Sequence operator+(char lhs, Sequence&& rhs) { return lhs + rhs.s; }
inline Sequence operator+(Sequence&& lhs, char rhs) { return lhs.s + rhs; }

inline bool operator==(const Sequence& lhs, const Sequence& rhs) { return lhs.s == rhs.s; }
inline bool operator==(const std::string& lhs, const Sequence& rhs) { return lhs == rhs.s; }
inline bool operator==(const Sequence& lhs, const std::string& rhs) { return lhs.s == rhs; }
inline bool operator==(const char* lhs, const Sequence& rhs) { return lhs == rhs.s; }
inline bool operator==(const Sequence& lhs, const char* rhs) { return lhs.s == rhs; }

inline bool operator!=(const Sequence& lhs, const Sequence& rhs) { return lhs.s != rhs.s; }
inline bool operator!=(const std::string& lhs, const Sequence& rhs) { return lhs != rhs.s; }
inline bool operator!=(const Sequence& lhs, const std::string& rhs) { return lhs.s != rhs; }
inline bool operator!=(const char* lhs, const Sequence& rhs) { return lhs != rhs.s; }
inline bool operator!=(const Sequence& lhs, const char* rhs) { return lhs.s != rhs; }

inline bool operator<(const Sequence& lhs, const Sequence& rhs) { return lhs.s < rhs.s; }
inline bool operator<(const std::string& lhs, const Sequence& rhs) { return lhs < rhs.s; }
inline bool operator<(const Sequence& lhs, const std::string& rhs) { return lhs.s < rhs; }
inline bool operator<(const char* lhs, const Sequence& rhs) { return lhs < rhs.s; }
inline bool operator<(const Sequence& lhs, const char* rhs) { return lhs.s < rhs; }

inline bool operator<=(const Sequence& lhs, const Sequence& rhs) { return lhs.s <= rhs.s; }
inline bool operator<=(const std::string& lhs, const Sequence& rhs) { return lhs <= rhs.s; }
inline bool operator<=(const Sequence& lhs, const std::string& rhs) { return lhs.s <= rhs; }
inline bool operator<=(const char* lhs, const Sequence& rhs) { return lhs <= rhs.s; }
inline bool operator<=(const Sequence& lhs, const char* rhs) { return lhs.s <= rhs; }

inline bool operator>(const Sequence& lhs, const Sequence& rhs) { return lhs.s > rhs.s; }
inline bool operator>(const std::string& lhs, const Sequence& rhs) { return lhs > rhs.s; }
inline bool operator>(const Sequence& lhs, const std::string& rhs) { return lhs.s > rhs; }
inline bool operator>(const char* lhs, const Sequence& rhs) { return lhs > rhs.s; }
inline bool operator>(const Sequence& lhs, const char* rhs) { return lhs.s > rhs; }

inline bool operator>=(const Sequence& lhs, const Sequence& rhs) { return lhs.s >= rhs.s; }
inline bool operator>=(const std::string& lhs, const Sequence& rhs) { return lhs >= rhs.s; }
inline bool operator>=(const Sequence& lhs, const std::string& rhs) { return lhs.s >= rhs; }
inline bool operator>=(const char* lhs, const Sequence& rhs) { return lhs >= rhs.s; }
inline bool operator>=(const Sequence& lhs, const char* rhs) { return lhs.s >= rhs; }

inline void swap(Sequence& x, Sequence& y) { swap(x.s, y.s); }

inline std::istream& operator>>(std::istream& is, Sequence& rhs) { is >> rhs.s; return is; }
inline std::ostream& operator<<(std::ostream& os, const Sequence& rhs) { os << rhs.s; return os; }

inline std::istream& getline(std::istream& is, Sequence& seq, char delim) { auto& ret = getline(is, seq.s, delim); seq.validate(); seq.capitalize(); return ret; }
inline std::istream& getline(std::istream&& is, Sequence& seq, char delim) { auto& ret = getline(is, seq.s, delim); seq.validate(); seq.capitalize(); return ret; }
inline std::istream& getline(std::istream& is, Sequence& seq) { auto& ret = getline(is, seq.s); seq.validate(); seq.capitalize(); return ret; }
inline std::istream& getline(std::istream&& is, Sequence& seq) { auto& ret = getline(is, seq.s); seq.validate(); seq.capitalize(); return ret; }

inline void Sequence::reverseComplement() {
	std::reverse(s.begin(), s.end());
	std::transform(s.begin(), s.end(), s.begin(),
        [] (char base) { return ~Base(base); }
    );
}

inline Sequence Sequence::getReverseComplement() {
    Sequence seq;
    seq.reserve(s.size());
    std::for_each(s.rbegin(), s.rend(),
        [&] (char base) { seq.s += ~Base(base); }
    );
    return seq;
}

inline Sequence Sequence::operator~() {
    return getReverseComplement();
}

inline void Sequence::validate(const std::string& bases) {
    std::for_each(bases.begin(), bases.end(),
        [] (char base) { Base::validate(base); }
    );
}

inline void Sequence::validate(const char* bases) {
    for (unsigned i = 0; bases[i] != '\0'; ++i) {
        Base::validate(bases[i]);  
    }
}

inline void Sequence::validate(std::initializer_list<char> bases) {
    std::for_each(bases.begin(), bases.end(),
        [] (char base) { Base::validate(base); }
    );
}

inline void Sequence::validate() { validate(s); }

inline void Sequence::capitalize() {
    std::for_each(s.begin(), s.end(),
        [] (char& base) { base = Base::capitalize(base); }
    );
}

inline void Sequence::set(std::string bases) { s = std::move(bases); }
inline void Sequence::set(const char* bases) { s = std::string(bases); }

} // namespace btl

#endif
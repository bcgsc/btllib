#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "basic_sequence.h"

namespace btl {

class Sequence: public BasicSequence {

public:

    using Base = BasicSequence::Base;

    Sequence() {}
    Sequence(const BasicSequence& seq): BasicSequence(seq) {}
    Sequence(const Sequence& seq): BasicSequence(seq) {}
    Sequence(std::string bases): BasicSequence(bases) {}
    Sequence(const Sequence& seq, size_t pos, size_t len = npos): BasicSequence(seq, pos, len) {}
    Sequence(const std::string& bases, size_t pos, size_t len = npos): BasicSequence(bases, pos, len) {}
    Sequence(const char* bases): BasicSequence(bases) {}
    Sequence(const char* bases, size_t n): BasicSequence(bases, n) {}
    Sequence(size_t n, char base): BasicSequence(n, base) {}
    template <class InputIterator>
    Sequence(InputIterator first, InputIterator last): BasicSequence(first, last) {}
    Sequence(std::initializer_list<char> bases): BasicSequence(bases) {}
    Sequence(BasicSequence&& seq): BasicSequence(seq) {}
    Sequence(Sequence&& seq): BasicSequence(seq) {}

    Sequence& operator=(const Sequence& seq) { BasicSequence::operator=(seq); return *this; }
    Sequence& operator=(std::string bases) { BasicSequence::operator=(bases); return *this; }
    Sequence& operator=(const char* bases) { BasicSequence::operator=(bases); return *this; }
    Sequence& operator=(std::initializer_list<char> bases) { BasicSequence::operator=(bases); return *this; }
    Sequence& operator=(BasicSequence&& seq) { BasicSequence::operator=(seq); return *this; }
    Sequence& operator=(Sequence&& seq) { BasicSequence::operator=(seq); return *this; }
    Sequence& operator=(std::string&& bases) { BasicSequence::operator=(bases); return *this; }

    void reverseComplement();
    Sequence getReverseComplement();
    Sequence operator~();

};

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

}

#endif
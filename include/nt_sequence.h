#ifndef NT_SEQUENCE_H
#define NT_SEQUENCE_H

#include "basic_nt_sequence.h"

namespace btl {

class NtSequence: public BasicNtSequence {

public:

    using Base = BasicNtSequence::Base;

    NtSequence() {}
    NtSequence(const BasicNtSequence& seq): BasicNtSequence(seq) {}
    NtSequence(const NtSequence& seq): BasicNtSequence(seq) {}
    NtSequence(std::string bases): BasicNtSequence(bases) {}
    NtSequence(const NtSequence& seq, size_t pos, size_t len = npos): BasicNtSequence(seq, pos, len) {}
    NtSequence(const std::string& bases, size_t pos, size_t len = npos): BasicNtSequence(bases, pos, len) {}
    NtSequence(const char* bases): BasicNtSequence(bases) {}
    NtSequence(const char* bases, size_t n): BasicNtSequence(bases, n) {}
    NtSequence(size_t n, char base): BasicNtSequence(n, base) {}
    template <class InputIterator>
    NtSequence(InputIterator first, InputIterator last): BasicNtSequence(first, last) {}
    NtSequence(std::initializer_list<char> bases): BasicNtSequence(bases) {}
    NtSequence(BasicNtSequence&& seq): BasicNtSequence(seq) {}
    NtSequence(NtSequence&& seq): BasicNtSequence(seq) {}

    NtSequence& operator=(const NtSequence& seq) { BasicNtSequence::operator=(seq); return *this; }
    NtSequence& operator=(std::string bases) { BasicNtSequence::operator=(bases); return *this; }
    NtSequence& operator=(const char* bases) { BasicNtSequence::operator=(bases); return *this; }
    NtSequence& operator=(std::initializer_list<char> bases) { BasicNtSequence::operator=(bases); return *this; }
    NtSequence& operator=(BasicNtSequence&& seq) { BasicNtSequence::operator=(seq); return *this; }
    NtSequence& operator=(NtSequence&& seq) { BasicNtSequence::operator=(seq); return *this; }
    NtSequence& operator=(std::string&& bases) { BasicNtSequence::operator=(bases); return *this; }

    void reverseComplement();
    NtSequence getReverseComplement();
    NtSequence operator~();

};

inline void NtSequence::reverseComplement() {
	std::reverse(s.begin(), s.end());
	std::transform(s.begin(), s.end(), s.begin(),
        [] (char base) { return ~Base(base); }
    );
}

inline NtSequence NtSequence::getReverseComplement() {
    NtSequence seq;
    seq.reserve(s.size());
    std::for_each(s.rbegin(), s.rend(),
        [&] (char base) { seq.s += ~Base(base); }
    );
    return seq;
}

inline NtSequence NtSequence::operator~() {
    return getReverseComplement();
}

}

#endif
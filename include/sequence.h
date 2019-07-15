#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "basic_sequence.h"

namespace btl {

class Sequence: public BasicSequence {

public:

    using Base = BasicSequence::Base;

    Sequence(const BasicSequence& seq): BasicSequence(seq) {}
    Sequence(const Sequence& seq): BasicSequence(seq) {}
    Sequence(BasicSequence&& seq): BasicSequence(seq) {}
    Sequence(Sequence&& seq): BasicSequence(seq) {}
    using BasicSequence::BasicSequence;

    Sequence& operator=(const BasicSequence& seq) { BasicSequence::operator=(seq); return *this; }
    Sequence& operator=(const Sequence& seq) { BasicSequence::operator=(seq); return *this; }
    Sequence& operator=(BasicSequence&& seq) { BasicSequence::operator=(seq); return *this; }
    Sequence& operator=(Sequence&& seq) { BasicSequence::operator=(seq); return *this; }
    using BasicSequence::operator=;

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
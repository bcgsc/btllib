#ifndef CONST_SEQUENCE_H
#define CONST_SEQUENCE_H

#include "basic_sequence.h"
#include "sequence.h"

namespace btl {

class ConstSequence: public BasicSequence {

public:

    using Base = BasicSequence::Base;

    ConstSequence(const BasicSequence& seq): BasicSequence(seq) {}
    ConstSequence(const ConstSequence& seq): BasicSequence(seq) {}
    ConstSequence(BasicSequence&& seq): BasicSequence(seq) {}
    ConstSequence(Sequence&& seq): BasicSequence(seq) {}
    ConstSequence(ConstSequence&& seq): BasicSequence(seq) {}
    using BasicSequence::BasicSequence;

};

}

#endif
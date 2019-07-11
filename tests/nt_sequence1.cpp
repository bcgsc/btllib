#include "nt_sequence.h"

#include <string>
#include <cstring>

int main() {
    btl::NtSequence seq1("BAR");
    assert(seq1 == std::string("BAR"));
    seq1 = "GATTACA";
    btl::NtSequence seq2 = std::string("CAT");
    btl::NtSequence seq3 = seq1 + seq2;
    assert(seq3 == "GATTACACAT");

    seq3 += "ATTACC";
    assert(seq3 == "GATTACACATATTACC");

    assert(~seq3 == "GGTAATATGTGTAATC");
    assert(~~seq3 == seq3);

    assert(seq3.size() == 16);
    assert(std::strlen(seq3) == 16);

    seq2 += 'H';
    seq2 += 'A';
    seq2 += 'T';
    assert(seq2 == btl::NtSequence("CATHAT"));

    return 0;
}
#include "sequence.h"

#include <string>
#include <cstring>

int main() {
    btl::Sequence seq1("BAR");
    assert(seq1 == std::string("BAR"));
    seq1 = "GATTACA";

    btl::Sequence seq2 = std::string("CAT");
    
    btl::Sequence seq3 = seq1 + seq2;
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
    assert(seq2 == btl::Sequence("CATHAT"));

    seq1 = "actg";
    seq2 = "ACTG";
    assert(seq1 == seq2);

    seq3 = "ac" + btl::Sequence("TG");
    assert(seq1 == seq3);

    return 0;
}
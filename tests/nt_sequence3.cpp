#include "nt_sequence.h"

int main() {
    btl::NtSequence seq = "ACTG";
    seq[2] = '8';
    return 0;
}
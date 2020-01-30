#include "../include/btllib/kmer_set.hpp"

#include <cassert>
#include <iostream>
#include <string>

int
main()
{
  std::string seq = "CACTATCGACGATCATTCGAGCATCAGCGACTG";
  std::string seq2 = "GTAGTACGATCAGCGACTATCGAGCTACGAGCA";
  assert(seq.size() == seq2.size());

  btllib::KmerSet kmer_set(seq.size() / 2, 1024 * 1024);
  kmer_set.insert(seq);
  assert(kmer_set.contains(seq) == (seq.size() - seq.size() / 2 + 1));
  assert(kmer_set.contains(seq2) <= 1);

  return 0;
}
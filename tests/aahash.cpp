#include "btllib/aahash.hpp"

#include "helpers.hpp"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>

int
main()
{
  const std::string seq = "RITMLYTIRITMLYTI";
  const unsigned k = 8, h = 3;
  btllib::AAHash aahash(seq, h, k);
  unsigned num_kmers = seq.size() - k + 2;
  aahash.roll();
  --num_kmers;
  uint64_t first_hash = aahash.hashes()[0];
  while (--num_kmers) {
    TEST_ASSERT(aahash.roll());
  }
  TEST_ASSERT(!aahash.roll());
  uint64_t last_hash = aahash.hashes()[0];
  TEST_ASSERT_EQ(first_hash, last_hash)
  return 0;
}
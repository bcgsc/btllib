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
  PRINT_TEST_NAME("Level 1 hashing")
  btllib::AAHash aahash(seq, h, k, 1);
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

  PRINT_TEST_NAME("Level 2 hashing")
  btllib::AAHash aahash2(seq, h, k, 2);
  num_kmers = seq.size() - k + 2;
  aahash2.roll();
  --num_kmers;
  first_hash = aahash2.hashes()[0];
  while (--num_kmers) {
    TEST_ASSERT(aahash2.roll());
  }
  TEST_ASSERT(!aahash2.roll());
  last_hash = aahash2.hashes()[0];
  TEST_ASSERT_EQ(first_hash, last_hash)

  PRINT_TEST_NAME("Level 3 hashing")
  btllib::AAHash aahash3(seq, h, k, 3);
  num_kmers = seq.size() - k + 2;
  aahash3.roll();
  --num_kmers;
  first_hash = aahash3.hashes()[0];
  while (--num_kmers) {
    TEST_ASSERT(aahash3.roll());
  }
  TEST_ASSERT(!aahash3.roll());
  last_hash = aahash3.hashes()[0];
  TEST_ASSERT_EQ(first_hash, last_hash)
  return 0;
}
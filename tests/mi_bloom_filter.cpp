#include "btllib/mi_bloom_filter.hpp"

#include <sdsl/bit_vector_il.hpp>
#include <sdsl/rank_support.hpp>

typedef uint16_t ID;

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <mutex>
#include <string>

int
main()
{
  std::cerr << "Testing MIBloomFilter" << std::endl;
  unsigned m_k = 10; // kmer length
  unsigned m_h = 4;  // # of hash function
  sdsl::bit_vector m_bv;
  std::vector<ID> m_counts;
  const std::vector<std::string>& m_spacedSeeds = { "10001", "11111" };

  btllib::MIBloomFilter<ID> mi_bf(m_h, m_k, m_bv, m_spacedSeeds);

  m_counts = std::vector<ID>(mi_bf.get_pop(), 0);
  return 0;
}

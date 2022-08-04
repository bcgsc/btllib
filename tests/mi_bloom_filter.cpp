#include "btllib/mi_bloom_filter.hpp"

#include "helpers.hpp"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>

int
main()
{
  std::cerr << "Testing multi-indexed BloomFilter" << std::endl;
  btllib::MIBloomFilter<uint8_t> mi_bf(1024 * 1024, 3, "ntHash");
  mi_bf.insert_bv({ 1, 10, 100 });
  mi_bf.insert_bv({ 100, 200, 300 });

  mi_bf.complete_bv_insertion();

  TEST_ASSERT(mi_bf.bv_contains({ 1, 10, 100 }));
  TEST_ASSERT(mi_bf.bv_contains({ 100, 200, 300 }));
  TEST_ASSERT(!mi_bf.bv_contains({ 1, 20, 100 }));

  uint8_t ID = 12;
  mi_bf.insert_id({ 1, 10, 100 }, ID);
  
  return 0;
}

#include "../include/btllib/bloom_filter.hpp"

#include "helpers.hpp"

#include <cstdio>

int
main()
{
  btllib::BloomFilter bf(1024 * 1024, 3);
  bf.insert({ 1, 10, 100 });
  bf.insert({ 100, 200, 300 });

  assert(bf.contains({ 1, 10, 100 }));
  assert(bf.contains({ 100, 200, 300 }));
  assert(!bf.contains({ 1, 20, 100 }));

  auto filename = get_random_name(64);
  bf.write(filename);

  btllib::BloomFilter bf2(filename);

  assert(bf2.contains({ 1, 10, 100 }));
  assert(bf2.contains({ 100, 200, 300 }));
  assert(!bf2.contains({ 1, 20, 100 }));

  std::remove(filename.c_str());

  return 0;
}
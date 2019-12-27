#ifndef BTLLIB_COUNTING_KMER_SET_HPP
#define BTLLIB_COUNTING_KMER_SET_HPP

#include "counting_bloom_filter.hpp"

#include <string>

namespace btllib {

class CountingKmerSet
{

public:
  CountingKmerSet(unsigned k, size_t bytes);

  void insert(const std::string& kmer);
  void insert(const char* kmer);

  unsigned count(const std::string& kmer);
  unsigned count(const char* kmer);

private:
  unsigned k;
  CountingBloomFilter bf;
};

inline CountingKmerSet::CountingKmerSet(unsigned k, size_t bytes)
  : k(k)
  , bf(bytes)
{}

inline void
CountingKmerSet::insert(const std::string& kmer)
{
  insert(kmer.c_str());
}

inline void
CountingKmerSet::insert(const char* kmer)
{}

inline unsigned
CountingKmerSet::count(const std::string& kmer)
{
  return count(kmer.c_str());
}

inline unsigned
CountingKmerSet::count(const char* kmer)
{
  return 0;
}

} // namespace btllib

#endif
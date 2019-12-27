#ifndef BTLLIB_COUNTING_KMER_SET_HPP
#define BTLLIB_COUNTING_KMER_SET_HPP

#include "counting_bloom_filter.hpp"

#include <string>

namespace btllib {

class CountingKmerSet
{

public:
  CountingKmerSet(unsigned k, size_t bytes);

  void insert(const std::string& seq);
  void insert(const char* seq);

  unsigned count(const std::string& seq);
  unsigned count(const char* seq);

private:
  unsigned k;
  CountingBloomFilter bf;
};

inline CountingKmerSet::CountingKmerSet(unsigned k, size_t bytes)
  : k(k)
  , bf(bytes)
{}

inline void
CountingKmerSet::insert(const std::string& seq)
{
  insert(seq.c_str());
}

inline void
CountingKmerSet::insert(const char* seq)
{}

inline unsigned
CountingKmerSet::count(const std::string& seq)
{
  return count(seq.c_str());
}

inline unsigned
CountingKmerSet::count(const char* seq)
{
  return 0;
}

} // namespace btllib

#endif
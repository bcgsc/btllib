#ifndef BTLLIB_KMER_SET_HPP
#define BTLLIB_KMER_SET_HPP

#include "bloom_filter.hpp"

#include <string>

namespace btllib {

class KmerSet
{

public:
  KmerSet(unsigned k, size_t bytes);

  void insert(const std::string& kmer);
  void insert(const char* kmer);

  bool contains(const std::string& kmer);
  bool contains(const char* kmer);

private:
  unsigned k;
  BloomFilter bf;
};

inline KmerSet::KmerSet(unsigned k, size_t bytes)
  : k(k)
  , bf(bytes)
{}

inline void
KmerSet::insert(const std::string& kmer)
{
  insert(kmer.c_str());
}

inline void
KmerSet::insert(const char* kmer)
{}

inline bool
KmerSet::contains(const std::string& kmer)
{
  return contains(kmer.c_str());
}
inline bool
KmerSet::contains(const char* kmer)
{
  return false;
}

} // namespace btllib

#endif
#ifndef BTLLIB_KMER_SET_HPP
#define BTLLIB_KMER_SET_HPP

#include "bloom_filter.hpp"

#include <string>

namespace btllib {

class KmerSet
{

public:
  KmerSet(unsigned k, size_t bytes);

  void insert(const std::string& seq);
  void insert(const char* seq);

  bool contains(const std::string& seq);
  bool contains(const char* seq);

private:
  unsigned k;
  BloomFilter bf;
  static const size_t HASH_NUM = 4;
};

inline KmerSet::KmerSet(unsigned k, size_t bytes)
  : k(k)
  , bf(bytes, HASH_NUM)
{}

inline void
KmerSet::insert(const std::string& seq)
{
  insert(seq.c_str());
}

inline void
KmerSet::insert(const char* seq)
{
  (void)seq;
}

inline bool
KmerSet::contains(const std::string& seq)
{
  return contains(seq.c_str());
}
inline bool
KmerSet::contains(const char* seq)
{
  (void)seq;
  return false;
}

} // namespace btllib

#endif
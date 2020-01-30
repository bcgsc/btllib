#ifndef BTLLIB_KMER_SET_HPP
#define BTLLIB_KMER_SET_HPP

#include "bloom_filter.hpp"
#include "rolling_hash.hpp"

#include <string>

namespace btllib {

class KmerSet
{

public:
  KmerSet(unsigned k, size_t bytes, unsigned hash_num = 4);

  void insert(const std::string& seq);
  void insert(const char* seq, size_t seq_len);

  unsigned contains(const std::string& seq);
  unsigned contains(const char* seq, size_t seq_len);

private:
  unsigned k;
  BloomFilter bf;
};

inline KmerSet::KmerSet(unsigned k, size_t bytes, unsigned hash_num)
  : k(k)
  , bf(bytes, hash_num)
{}

inline void
KmerSet::insert(const std::string& seq)
{
  insert(seq.c_str(), seq.size());
}

inline void
KmerSet::insert(const char* seq, size_t seq_len)
{
  RollingHash rolling_hash(seq, seq_len, k, bf.get_hash_num());
  while (rolling_hash.roll()) {
    bf.insert(rolling_hash.hashes());
  }
}

inline unsigned
KmerSet::contains(const std::string& seq)
{
  return contains(seq.c_str(), seq.size());
}
inline unsigned
KmerSet::contains(const char* seq, size_t seq_len)
{
  unsigned count = 0;
  RollingHash rolling_hash(seq, seq_len, k, bf.get_hash_num());
  while (rolling_hash.roll()) {
    if (bf.contains(rolling_hash.hashes())) {
      count++;
    }
  }
  return count;
}

} // namespace btllib

#endif
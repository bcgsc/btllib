#ifndef BTLLIB_COUNTING_BLOOM_FILTER_HPP
#define BTLLIB_COUNTING_BLOOM_FILTER_HPP

#include <cstdint>
#include <vector>

namespace btllib {

class CountingBloomFilter
{

public:
  CountingBloomFilter(size_t size);

  void insert(const std::vector<uint64_t>& hashes);
  void insert(const uint64_t* hashes);

  unsigned count(const std::vector<uint64_t>& hashes);
  unsigned count(const uint64_t* hashes);

private:
  size_t size;
};

inline CountingBloomFilter::CountingBloomFilter(size_t size)
  : size(size)
{}

inline void
CountingBloomFilter::insert(const std::vector<uint64_t>& hashes)
{
  insert(hashes.data());
}

inline void
CountingBloomFilter::insert(const uint64_t* hashes)
{
  (void)hashes;
}

inline unsigned
CountingBloomFilter::count(const std::vector<uint64_t>& hashes)
{
  return count(hashes.data());
}
inline unsigned
CountingBloomFilter::count(const uint64_t* hashes)
{
  (void)hashes;
  return 0;
}

} // namespace btllib

#endif
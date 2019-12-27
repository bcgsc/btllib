#ifndef BTLLIB_BLOOM_FILTER_HPP
#define BTLLIB_BLOOM_FILTER_HPP

#include <cstdint>
#include <vector>

namespace btllib {

class BloomFilter
{

public:
  BloomFilter(size_t size);

  void insert(const std::vector<uint64_t>& hashes);
  void insert(const uint64_t* hashes);

  bool contains(const std::vector<uint64_t>& hashes);
  bool contains(const uint64_t* hashes);

private:
  size_t size;
};

inline BloomFilter::BloomFilter(size_t size): size(size) {}

inline void
BloomFilter::insert(const std::vector<uint64_t>& hashes)
{}

inline void
BloomFilter::insert(const uint64_t* hashes)
{}

inline bool
BloomFilter::contains(const std::vector<uint64_t>& hashes)
{
    return contains(hashes.data());
}
inline bool
BloomFilter::contains(const uint64_t* hashes)
{
    return false;
}

} // namespace btllib

#endif
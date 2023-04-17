#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"
#include "btllib/nthash.hpp"
#include "btllib/status.hpp"

#include "cpptoml.h"

#include <atomic>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace btllib {

  using CountingBloomFilter8 = CountingBloomFilter<uint8_t>;
  using CountingBloomFilter16 = CountingBloomFilter<uint16_t>;
  using CountingBloomFilter32 = CountingBloomFilter<uint32_t>;

  using KmerCountingBloomFilter8 = KmerCountingBloomFilter<uint8_t>;
  using KmerCountingBloomFilter16 = KmerCountingBloomFilter<uint16_t>;
  using KmerCountingBloomFilter32 = KmerCountingBloomFilter<uint32_t>;

  template<typename T>
  inline CountingBloomFilter<T>::CountingBloomFilter(size_t bytes,
              unsigned hash_num,
              std::string hash_fn)
  : bytes(
  size_t(std::ceil(double(bytes) / sizeof(uint64_t)) * sizeof(uint64_t)))
  , array_size(get_bytes() / sizeof(array[0]))
  , hash_num(hash_num)
  , hash_fn(std::move(hash_fn))
  , array(new std::atomic<T>[array_size])
  {
    check_error(bytes == 0, "CountingBloomFilter: memory budget must be >0!");
    check_error(hash_num == 0,
      "CountingBloomFilter: number of hash values must be >0!");
    check_error(
    hash_num > MAX_HASH_VALUES,
    "CountingBloomFilter: number of hash values cannot be over 1024!");
    check_warning(sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
        "Atomic primitives take extra memory. CountingBloomFilter will "
        "have less than " +
        std::to_string(bytes) + " for bit array.");
    std::memset((void*)array.get(), 0, array_size * sizeof(array[0]));
  }

  /*
  * Assumes min_count is not std::numeric_limits<T>::max()
  */
  template<typename T>
  inline void
  CountingBloomFilter<T>::insert(const uint64_t* hashes, T min_val)
  {
    // Update flag to track if increment is done on at least one counter
    bool update_done = false;
    T new_val, tmp_min_val;
    while (true) {
    new_val = min_val + 1;
    for (size_t i = 0; i < hash_num; ++i) {
      tmp_min_val = min_val;
      update_done |= array[hashes[i] % array_size].compare_exchange_strong(
        tmp_min_val, new_val);
    }
    if (update_done) {
      break;
    }
    min_val = contains(hashes);
      if (min_val == std::numeric_limits<T>::max()) {
        break;
      }
    }
  }

  template<typename T>
  inline void
  CountingBloomFilter<T>::insert(const uint64_t* hashes)
  {
    contains_insert(hashes);
  }

  template<typename T>
  inline void
  CountingBloomFilter<T>::remove(const uint64_t* hashes)
  {
    // Update flag to track if increment is done on at least one counter
    bool update_done = false;
    T min_val = contains(hashes);
    T new_val, tmp_min_val;
    while (true) {
      new_val = min_val - 1;
      for (size_t i = 0; i < hash_num; ++i) {
        tmp_min_val = min_val;
        update_done |= array[hashes[i] % array_size].compare_exchange_strong(
          tmp_min_val, new_val);
      }
      if (update_done) {
        break;
      }
      min_val = contains(hashes);
      if (min_val == std::numeric_limits<T>::max()) {
        break;
      }
    }
  }

  template<typename T>
  inline void
  CountingBloomFilter<T>::clear(const uint64_t* hashes)
  {
    // Update flag to track if increment is done on at least one counter
    bool update_done = false;
    T min_val = contains(hashes);
    T new_val, tmp_min_val;
    while (true) {
      new_val = 0;
      for (size_t i = 0; i < hash_num; ++i) {
        tmp_min_val = min_val;
        update_done |= array[hashes[i] % array_size].compare_exchange_strong(
          tmp_min_val, new_val);
      }
      if (update_done) {
        break;
      }
      min_val = contains(hashes);
      if (min_val == std::numeric_limits<T>::max()) {
        break;
      }
    }
  }

  template<typename T>
  inline T
  CountingBloomFilter<T>::contains(const uint64_t* hashes) const
  {
    T min = array[hashes[0] % array_size];
    for (size_t i = 1; i < hash_num; ++i) {
      const size_t idx = hashes[i] % array_size;
      if (array[idx] < min) {
        min = array[idx];
      }
    }
    return min;
  }

  template<typename T>
  inline T
  CountingBloomFilter<T>::contains_insert(const uint64_t* hashes)
  {
    const auto count = contains(hashes);
    if (count < std::numeric_limits<T>::max()) {
      insert(hashes, count);
    }
    return count;
  }

  template<typename T>
  inline T
  CountingBloomFilter<T>::insert_contains(const uint64_t* hashes)
  {
  const auto count = contains(hashes);
  if (count < std::numeric_limits<T>::max()) {
  insert(hashes, count);
  return count + 1;
  }
  return std::numeric_limits<T>::max();
  }

  template<typename T>
  inline T
  CountingBloomFilter<T>::insert_thresh_contains(const uint64_t* hashes,
            const T threshold)
  {
  const auto count = contains(hashes);
  if (count < threshold) {
  insert(hashes, count);
  return count + 1;
  }
  return count;
  }

  template<typename T>
  inline T
  CountingBloomFilter<T>::contains_insert_thresh(const uint64_t* hashes,
            const T threshold)
  {
  const auto count = contains(hashes);
  if (count < threshold) {
  insert(hashes, count);
  }
  return count;
  }

  template<typename T>
  inline uint64_t
  CountingBloomFilter<T>::get_pop_cnt(const T threshold) const
  {
  uint64_t pop_cnt = 0;
  // OpenMP make up your mind man. Using default(none) here causes errors on
  // some compilers and not others.
  // NOLINTNEXTLINE(openmp-use-default-none,-warnings-as-errors)
  #pragma omp parallel for reduction(+ : pop_cnt)
  for (size_t i = 0; i < array_size; ++i) {
  if (array[i] >= threshold) {
  ++pop_cnt;
  }
  }
  return pop_cnt;
  }

  template<typename T>
  inline double
  CountingBloomFilter<T>::get_occupancy(const T threshold) const
  {
  return double(get_pop_cnt(threshold)) / double(array_size);
  }

  template<typename T>
  inline double
  CountingBloomFilter<T>::get_fpr(const T threshold) const
  {
  return std::pow(get_occupancy(threshold), double(hash_num));
  }

  template<typename T>
  inline CountingBloomFilter<T>::CountingBloomFilter(const std::string& path)
  : CountingBloomFilter<T>::CountingBloomFilter(
  std::make_shared<BloomFilterInitializer>(path,
            COUNTING_BLOOM_FILTER_SIGNATURE))
  {
  }

  template<typename T>
  inline CountingBloomFilter<T>::CountingBloomFilter(
  const std::shared_ptr<BloomFilterInitializer>& bfi)
  : bytes(*bfi->table->get_as<decltype(bytes)>("bytes"))
  , array_size(bytes / sizeof(array[0]))
  , hash_num(*(bfi->table->get_as<decltype(hash_num)>("hash_num")))
  , hash_fn(bfi->table->contains("hash_fn")
    ? *(bfi->table->get_as<decltype(hash_fn)>("hash_fn"))
    : "")
  , array(new std::atomic<T>[array_size])
  {
  check_warning(sizeof(uint8_t) != sizeof(std::atomic<uint8_t>),
      "Atomic primitives take extra memory. CountingBloomFilter will "
      "have less than " +
      std::to_string(bytes) + " for bit array.");
  const auto loaded_counter_bits =
  *(bfi->table->get_as<size_t>("counter_bits"));
  check_error(sizeof(array[0]) * CHAR_BIT != loaded_counter_bits,
    "CountingBloomFilter" +
      std::to_string(sizeof(array[0]) * CHAR_BIT) +
      " tried to load a file of CountingBloomFilter" +
      std::to_string(loaded_counter_bits));
  bfi->ifs.read((char*)array.get(),
      std::streamsize(array_size * sizeof(array[0])));
  }

  template<typename T>
  inline void
  CountingBloomFilter<T>::save(const std::string& path)
  {
  /* Initialize cpptoml root table
  Note: Tables and fields are unordered
  Ordering of table is maintained by directing the table
  to the output stream immediately after completion  */
  auto root = cpptoml::make_table();

  /* Initialize bloom filter section and insert fields
  and output to ostream */
  auto header = cpptoml::make_table();
  header->insert("bytes", get_bytes());
  header->insert("hash_num", get_hash_num());
  if (!hash_fn.empty()) {
  header->insert("hash_fn", hash_fn);
  }
  header->insert("counter_bits", size_t(sizeof(array[0]) * CHAR_BIT));
  std::string header_string = COUNTING_BLOOM_FILTER_SIGNATURE;
  header_string =
  header_string.substr(1, header_string.size() - 2); // Remove [ ]
  root->insert(header_string, header);

  BloomFilter::save(
  path, *root, (char*)array.get(), array_size * sizeof(array[0]));
  }

  template<typename T>
  inline KmerCountingBloomFilter<T>::KmerCountingBloomFilter(size_t bytes,
                unsigned hash_num,
                unsigned k)
  : k(k)
  , counting_bloom_filter(bytes, hash_num, HASH_FN)
  {
  }

  template<typename T>
  inline void
  KmerCountingBloomFilter<T>::insert(const char* seq, size_t seq_len)
  {
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
  counting_bloom_filter.insert(nthash.hashes());
  }
  }

  template<typename T>
  inline void
  KmerCountingBloomFilter<T>::remove(const char* seq, size_t seq_len)
  {
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
  counting_bloom_filter.remove(nthash.hashes());
  }
  }

  template<typename T>
  inline void
  KmerCountingBloomFilter<T>::clear(const char* seq, size_t seq_len)
  {
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
  counting_bloom_filter.clear(nthash.hashes());
  }
  }

  template<typename T>
  inline uint64_t
  KmerCountingBloomFilter<T>::contains(const char* seq, size_t seq_len) const
  {
  uint64_t sum = 0;
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
  sum += counting_bloom_filter.contains(nthash.hashes());
  }
  return sum;
  }

  template<typename T>
  inline T
  KmerCountingBloomFilter<T>::contains_insert(const char* seq, size_t seq_len)
  {
  uint64_t sum = 0;
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
  sum += counting_bloom_filter.contains_insert(nthash.hashes());
  }
  return sum;
  }

  template<typename T>
  inline T
  KmerCountingBloomFilter<T>::insert_contains(const char* seq, size_t seq_len)
  {
  uint64_t sum = 0;
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
  sum += counting_bloom_filter.insert_contains(nthash.hashes());
  }
  return sum;
  }

  template<typename T>
  inline T
  KmerCountingBloomFilter<T>::insert_thresh_contains(const char* seq,
              size_t seq_len,
              const T threshold)
  {
  uint64_t sum = 0;
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
  sum +=
  counting_bloom_filter.insert_thresh_contains(nthash.hashes(), threshold);
  }
  return sum;
  }

  template<typename T>
  inline T
  KmerCountingBloomFilter<T>::contains_insert_thresh(const char* seq,
              size_t seq_len,
              const T threshold)
  {
  uint64_t sum = 0;
  NtHash nthash(seq, seq_len, get_hash_num(), get_k());
  while (nthash.roll()) {
  sum +=
  counting_bloom_filter.contains_insert_thresh(nthash.hashes(), threshold);
  }
  return sum;
  }

  template<typename T>
  inline KmerCountingBloomFilter<T>::KmerCountingBloomFilter(
  const std::string& path)
  : KmerCountingBloomFilter<T>::KmerCountingBloomFilter(
  std::make_shared<BloomFilterInitializer>(
    path,
    KMER_COUNTING_BLOOM_FILTER_SIGNATURE))
  {
  }

  template<typename T>
  inline KmerCountingBloomFilter<T>::KmerCountingBloomFilter(
  const std::shared_ptr<BloomFilterInitializer>& bfi)
  : k(*(bfi->table->get_as<decltype(k)>("k")))
  , counting_bloom_filter(bfi)
  {
  check_error(counting_bloom_filter.hash_fn != HASH_FN,
    "KmerCountingBloomFilter: loaded hash function (" +
      counting_bloom_filter.hash_fn +
      ") is different from the one used by default (" + HASH_FN +
      ").");
  }

  template<typename T>
  inline void
  KmerCountingBloomFilter<T>::save(const std::string& path)
  {
  /* Initialize cpptoml root table
  Note: Tables and fields are unordered
  Ordering of table is maintained by directing the table
  to the output stream immediately after completion  */
  auto root = cpptoml::make_table();

  /* Initialize bloom filter section and insert fields
  and output to ostream */
  auto header = cpptoml::make_table();
  header->insert("bytes", get_bytes());
  header->insert("hash_num", get_hash_num());
  header->insert("hash_fn", get_hash_fn());
  header->insert("counter_bits",
      size_t(sizeof(counting_bloom_filter.array[0]) * CHAR_BIT));
  header->insert("k", k);
  std::string header_string = KMER_COUNTING_BLOOM_FILTER_SIGNATURE;
  header_string =
  header_string.substr(1, header_string.size() - 2); // Remove [ ]
  root->insert(header_string, header);

  BloomFilter::save(path,
      *root,
      (char*)counting_bloom_filter.array.get(),
      counting_bloom_filter.array_size *
      sizeof(counting_bloom_filter.array[0]));
  }
} // namespace btllib
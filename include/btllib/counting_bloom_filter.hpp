#ifndef BTLLIB_COUNTING_BLOOM_FILTER_HPP
#define BTLLIB_COUNTING_BLOOM_FILTER_HPP

#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"
#include "btllib/nthash.hpp"
#include "btllib/status.hpp"

// clang-format off
#include <limits>
#include "cpptoml.h"
// clang-format on

#include <atomic>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace btllib {

// NOLINTNEXTLINE(clang-diagnostic-unneeded-internal-declaration)
static const char* const COUNTING_BLOOM_FILTER_SIGNATURE =
  "[BTLCountingBloomFilter_v5]";
// NOLINTNEXTLINE(clang-diagnostic-unneeded-internal-declaration)
static const char* const KMER_COUNTING_BLOOM_FILTER_SIGNATURE =
  "[BTLKmerCountingBloomFilter_v5]";

template<typename T>
class KmerCountingBloomFilter;

/**
 * Counting Bloom filter data structure. Provides CountingBloomFilter8,
 * CountingBloomFilter16, and CountingBloomFilter32 classes with corresponding
 * bit-size counters.
 */
template<typename T>
class CountingBloomFilter
{

public:
  /** Construct a dummy k-mer Bloom filter (e.g. as a default argument). */
  CountingBloomFilter() {}

  /**
   * Construct an empty Counting Bloom filter of given size.
   *
   * @param bytes Filter size in bytes.
   * @param hash_num Number of hash values per element.
   * @param hash_fn Name of the hash function used. Used for metadata. Optional.
   */
  CountingBloomFilter(size_t bytes,
                      unsigned hash_num,
                      std::string hash_fn = "");

  /**
   * Load a Counting Bloom filter from a file.
   *
   * @param path Filepath to load from.
   */
  explicit CountingBloomFilter(const std::string& path);

  CountingBloomFilter(const CountingBloomFilter&) = delete;
  CountingBloomFilter(CountingBloomFilter&&) = delete;

  CountingBloomFilter& operator=(const CountingBloomFilter&) = delete;
  CountingBloomFilter& operator=(CountingBloomFilter&&) = delete;

  /**
   * Insert an element.
   *
   * @param hashes Integer array of the element's hash values. Array size should
   * equal the hash_num argument used when the Bloom filter was constructed.
   */
  void insert(const uint64_t* hashes);

  /**
   * Insert an element.
   *
   * @param hashes Integer vector of the element's hash values.
   */
  void insert(const std::vector<uint64_t>& hashes) { insert(hashes.data()); }

  /**
   * Delete an element.
   *
   * @param hashes Integer array of the element's hash values. Array size should
   * equal the hash_num argument used when the Bloom filter was constructed.
   */
  void remove(const uint64_t* hashes);

  /**
   * Delete an element.
   *
   * @param hashes Integer vector of the element's hash values.
   */
  void remove(const std::vector<uint64_t>& hashes) { remove(hashes.data()); }

  /**
   * Set the count of an element to zero.
   *
   * @param hashes Integer array of the element's hash values. Array size should
   * equal the hash_num argument used when the Bloom filter was constructed.
   */
  void clear(const uint64_t* hashes);

  /**
   * Set the count of an element to zero.
   *
   * @param hashes Integer vector of the element's hash values.
   */
  void clear(const std::vector<uint64_t>& hashes) { clear(hashes.data()); }

  /**
   * Get the count of an element.
   *
   * @param hashes Integer array of the element's hash values. Array size should
   * equal the hash_num argument used when the Bloom filter was constructed.
   *
   * @return The count of the queried element.
   */
  T contains(const uint64_t* hashes) const;

  /**
   * Get the count of an element.
   *
   * @param hashes Integer vector of the element's hash values.
   *
   * @return The count of the queried element.
   */
  T contains(const std::vector<uint64_t>& hashes) const
  {
    return contains(hashes.data());
  }

  /**
   * Get the count of an element and then increment the count.
   *
   * @param hashes Integer array of the element's hash values. Array size should
   * equal the hash_num argument used when the Bloom filter was constructed.
   *
   * @return The count of the queried element before insertion.
   */
  T contains_insert(const uint64_t* hashes);

  /**
   * Get the count of an element and then increment the count.
   *
   * @param hashes Integer vector of the element's hash values.
   *
   * @return The count of the queried element before insertion.
   */
  T contains_insert(const std::vector<uint64_t>& hashes)
  {
    return contains_insert(hashes.data());
  }

  /**
   * Increment an element's count and then return the count.
   *
   * @param hashes Integer array of the element's hash values. Array size
   * should equal the hash_num argument used when the Bloom filter was
   * constructed.
   *
   * @return The count of the queried element after insertion.
   */
  T insert_contains(const uint64_t* hashes);

  /**
   * Increment an element's count and then return the count.
   *
   * @param hashes Integer vector of the element's hash values.
   *
   * @return The count of the queried element after insertion.
   */
  T insert_contains(const std::vector<uint64_t>& hashes)
  {
    return insert_contains(hashes.data());
  }

  /**
   * Increment an element's count if it's not above the threshold and then
   * return the count.
   *
   * @param hashes Integer array of the element's hash values. Array size
   * should equal the hash_num argument used when the Bloom filter was
   * constructed.
   * @param threshold The threshold.
   *
   * @return The count of the queried element after insertion.
   */
  T insert_thresh_contains(const uint64_t* hashes, T threshold);

  /**
   * Increment an element's count if it's not above the threshold and then
   * return the count.
   *
   * @param hashes Integer array of the element's hash values. Array size
   * should equal the hash_num argument used when the Bloom filter was
   * constructed.
   * @param threshold The threshold.
   *
   * @return The count of the queried element after insertion.
   */
  T insert_thresh_contains(const std::vector<uint64_t>& hashes,
                           const T threshold)
  {
    return insert_thresh_contains(hashes.data(), threshold);
  }

  /**
   * Get the count of an element and then increment the count if it's not
   * above the threshold.
   *
   * @param hashes Integer array of the element's hash values. Array size
   * should equal the hash_num argument used when the Bloom filter was
   * constructed.
   * @param threshold The threshold.
   *
   * @return The count of the queried element before insertion.
   */
  T contains_insert_thresh(const uint64_t* hashes, T threshold);

  /**
   * Get the count of an element and then increment the count if it's not
   * above the threshold.
   *
   * @param hashes Integer vector of the element's hash values.
   * @param threshold The threshold.
   *
   * @return The count of the queried element before insertion.
   */
  T contains_insert_thresh(const std::vector<uint64_t>& hashes,
                           const T threshold)
  {
    return contains_insert_thresh(hashes.data(), threshold);
  }

  /** Get filter size in bytes. */
  size_t get_bytes() const { return bytes; }
  /** Get population count, i.e. the number of counters >= threshold in the
   * filter. */
  uint64_t get_pop_cnt(T threshold = 1) const;
  /** Get the fraction of the filter occupied by >= threshold counters. */
  double get_occupancy(T threshold = 1) const;
  /** Get the number of hash values per element. */
  unsigned get_hash_num() const { return hash_num; }
  /** Get the query false positive rate for elements with count >= threshold.
   *
   * @param threshold The threshold.
   *
   * @return The false positive rate.
   */
  double get_fpr(T threshold = 1) const;
  /** Get the name of the hash function used. */
  const std::string& get_hash_fn() const { return hash_fn; }

  /**
   * Save the Bloom filter to a file that can be loaded in the future.
   *
   * @param path Filepath to store filter at.
   */
  void save(const std::string& path);

  /**
   * Check whether the file at the given path is a saved Counting Bloom filter.
   *
   * @param path Filepath to check.
   */
  static bool is_bloom_file(const std::string& path)
  {
    return btllib::BloomFilter::check_file_signature(
      path, COUNTING_BLOOM_FILTER_SIGNATURE);
  }

private:
  CountingBloomFilter(const std::shared_ptr<BloomFilterInitializer>& bfi);

  void insert(const uint64_t* hashes, T min_val);

  friend class KmerCountingBloomFilter<T>;

  size_t bytes = 0;
  size_t array_size = 0;
  unsigned hash_num = 0;
  std::string hash_fn;
  std::unique_ptr<std::atomic<T>[]> array;
};

/**
 * Counting Bloom filter data structure that stores k-mers. Provides
 * KmerCountingBloomFilter8, KmerCountingBloomFilter16, and
 * KmerCountingBloomFilter32 classes with corresponding bit-size counters.
 */
template<typename T>
class KmerCountingBloomFilter
{

public:
  /** Construct a dummy k-mer Bloom filter (e.g. as a default argument). */
  KmerCountingBloomFilter() {}

  /**
   * Construct an empty k-mer Counting Bloom filter of given size.
   *
   * @param bytes Filter size in bytes.
   * @param hash_num Number of hash values per element.
   * @param k K-mer size.
   */
  KmerCountingBloomFilter(size_t bytes, unsigned hash_num, unsigned k);

  /**
   * Load a k-mer Counting Bloom filter from a file.
   *
   * @param path Filepath to load from.
   */
  explicit KmerCountingBloomFilter(const std::string& path);

  KmerCountingBloomFilter(const KmerCountingBloomFilter&) = delete;
  KmerCountingBloomFilter(KmerCountingBloomFilter&&) = delete;

  KmerCountingBloomFilter& operator=(const KmerCountingBloomFilter&) = delete;
  KmerCountingBloomFilter& operator=(KmerCountingBloomFilter&&) = delete;

  /**
   * Insert a sequence's k-mers into the filter.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   */
  void insert(const char* seq, size_t seq_len);

  /**
   * Insert a sequence's k-mers into the filter.
   *
   * @param seq Sequence to k-merize.
   */
  void insert(const std::string& seq) { insert(seq.c_str(), seq.size()); }

  /**
   * Insert a k-mer into the filter.
   *
   * @param hashes Integer array of the k-mer's hash values. Array size should
   * equal the hash_num argument used when the Bloom filter was constructed.
   */
  void insert(const uint64_t* hashes) { counting_bloom_filter.insert(hashes); }

  /**
   * Insert a k-mer into the filter.
   *
   * @param hashes Integer vector of the k-mer's hash values.
   */
  void insert(const std::vector<uint64_t>& hashes)
  {
    counting_bloom_filter.insert(hashes);
  }

  /**
   * Decrease the counts of a sequence's k-mers from the filter.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   */
  void remove(const char* seq, size_t seq_len);

  /**
   * Decrease the counts of a sequence's k-mers from the filter.
   *
   * @param seq Sequence to k-merize.
   */
  void remove(const std::string& seq) { remove(seq.c_str(), seq.size()); }

  /**
   * Decrease the counts of a sequence's k-mers from the filter.
   *
   * @param hashes Integer array of the k-mer's hash values. Array size should
   * equal the hash_num argument used when the Bloom filter was constructed.
   */
  void remove(const uint64_t* hashes) { counting_bloom_filter.remove(hashes); }

  /**
   * Decrease the counts of a sequence's k-mers from the filter.
   *
   * @param hashes Integer vector of the k-mer's hash values.
   */
  void remove(const std::vector<uint64_t>& hashes)
  {
    counting_bloom_filter.remove(hashes);
  }

  /**
   * Set the counts of a sequence's k-mers to zero in the filter.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   */
  void clear(const char* seq, size_t seq_len);

  /**
   * Set the counts of a sequence's k-mers to zero in the filter.
   *
   * @param seq Sequence to k-merize.
   */
  void clear(const std::string& seq) { clear(seq.c_str(), seq.size()); }

  /**
   * Set the counts of a sequence's k-mers to zero in the filter.
   *
   * @param hashes Integer array of the k-mer's hash values. Array size should
   * equal the hash_num argument used when the Bloom filter was constructed.
   */
  void clear(const uint64_t* hashes) { counting_bloom_filter.clear(hashes); }

  /**
   * Set the counts of a sequence's k-mers to zero in the filter.
   *
   * @param hashes Integer vector of the k-mer's hash values.
   */
  void clear(const std::vector<uint64_t>& hashes)
  {
    counting_bloom_filter.clear(hashes);
  }

  /**
   * Query the counts of k-mers of a sequence.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   *
   * @return The sum of counts of seq's k-mers found in the filter.
   */
  uint64_t contains(const char* seq, size_t seq_len) const;

  /**
   * Query the counts of k-mers of a sequence.
   *
   * @param seq Sequence to k-merize.
   *
   * @return The sum of counts of seq's k-mers found in the filter.
   */
  uint64_t contains(const std::string& seq) const
  {
    return contains(seq.c_str(), seq.size());
  }

  /**
   * Get a k-mer's count.
   *
   * @param hashes Integer array of k-mer's hash values. Array size should
   * equal the hash_num argument used when the Bloom filter was constructed.
   *
   * @return The count of the queried k-mer.
   */
  T contains(const uint64_t* hashes) const
  {
    return counting_bloom_filter.contains(hashes);
  }

  /**
   * Get a k-mer's count.
   *
   * @param hashes Integer vector of k-mer's hash values.
   *
   * @return The count of the queried k-mer.
   */
  T contains(const std::vector<uint64_t>& hashes) const
  {
    return counting_bloom_filter.contains(hashes);
  }

  /**
   * Get the counts of a sequence's k-mers and then increment the counts.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   *
   * @return The sum of counts of the queried k-mers before insertion.
   */
  T contains_insert(const char* seq, size_t seq_len);

  /**
   * Get the counts of a sequence's k-mers and then increment the counts.
   *
   * @param seq Sequence to k-merize.
   *
   * @return The sum of counts of the queried k-mers before insertion.
   */
  T contains_insert(const std::string& seq)
  {
    return contains_insert(seq.c_str(), seq.size());
  }

  /**
   * Get the count of a k-mer and then increment the count.
   *
   * @param hashes Integer array of the k-mers's hash values. Array size should
   * equal the hash_num argument used when the Bloom filter was constructed.
   *
   * @return The count of the queried k-mer before insertion.
   */
  T contains_insert(const uint64_t* hashes)
  {
    return counting_bloom_filter.contains_insert(hashes);
  }

  /**
   * Get the count of a k-mer and then increment the count.
   *
   * @param hashes Integer vector of the k-mer's hash values.
   *
   * @return The count of the queried k-mer before insertion.
   */
  T contains_insert(const std::vector<uint64_t>& hashes)
  {
    return counting_bloom_filter.contains_insert(hashes);
  }

  /**
   * Increment the counts of a sequence's k-mers and then return the counts.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   *
   * @return The sum of counts of the queried k-mers after insertion.
   */
  T insert_contains(const char* seq, size_t seq_len);

  /**
   * Increment the counts of a sequence's k-mers and then return the counts.
   *
   * @param seq Sequence to k-merize.
   *
   * @return The sum of counts of the queried k-mers after insertion.
   */
  T insert_contains(const std::string& seq)
  {
    return insert_contains(seq.c_str(), seq.size());
  }

  /**
   * Increment a k-mer's count and then return the count.
   *
   * @param hashes Integer array of the k-mer's hash values. Array size
   * should equal the hash_num argument used when the Bloom filter was
   * constructed.
   *
   * @return The count of the queried k-mer after insertion.
   */
  T insert_contains(const uint64_t* hashes)
  {
    return counting_bloom_filter.insert_contains(hashes);
  }

  /**
   * Increment a k-mer's count and then return the count.
   *
   * @param hashes Integer vector of the k-mer's hash values.
   *
   * @return The count of the queried k-mer after insertion.
   */
  T insert_contains(const std::vector<uint64_t>& hashes)
  {
    return counting_bloom_filter.insert_contains(hashes);
  }

  /**
   * Increment the counts of a sequence's k-mers if they are not above the
   * threshold and then return the counts.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   * @param threshold The threshold.
   *
   * @return The sum of counts of the queried k-mers after insertion.
   */
  T insert_thresh_contains(const char* seq, size_t seq_len, T threshold);

  /**
   * Increment the counts of a sequence's k-mers if they are not above the
   * threshold and then return the counts.
   *
   * @param seq Sequence to k-merize.
   * @param threshold The threshold.
   *
   * @return The sum of counts of the queried k-mers after insertion.
   */
  T insert_thresh_contains(const std::string& seq, const T threshold)
  {
    return insert_thresh_contains(seq.c_str(), seq.size(), threshold);
  }

  /**
   * Increment a k-mer's count if it's not above the threshold and then
   * return the count.
   *
   * @param hashes Integer array of the k-mer's hash values. Array size
   * should equal the hash_num argument used when the Bloom filter was
   * constructed.
   * @param threshold The threshold.
   *
   * @return The count of the queried k-mer after insertion.
   */
  T insert_thresh_contains(const uint64_t* hashes, const T threshold)
  {
    return counting_bloom_filter.insert_thresh_contains(hashes, threshold);
  }

  /**
   * Increment a k-mer's count if it's not above the threshold and then
   * return the count.
   *
   * @param hashes Integer array of the k-mer's hash values. Array size
   * should equal the hash_num argument used when the Bloom filter was
   * constructed.
   * @param threshold The threshold.
   *
   * @return The count of the queried k-mer after insertion.
   */
  T insert_thresh_contains(const std::vector<uint64_t>& hashes,
                           const T threshold)
  {
    return counting_bloom_filter.insert_thresh_contains(hashes, threshold);
  }

  /**
   * Get the counts of a sequence's k-mer's and then increment the counts if
   * they are not above the threshold.
   *
   * @param seq Sequence to k-merize.
   * @param seq_len Length of seq.
   * @param threshold The threshold.
   *
   * @return The sum of counts of the queried k-mers before insertion.
   */
  T contains_insert_thresh(const char* seq, size_t seq_len, T threshold);

  /**
   * Get the counts of a sequence's k-mer's and then increment the counts if
   * they are not above the threshold.
   *
   * @param seq Sequence to k-merize.
   * @param threshold The threshold.
   *
   * @return The sum of counts of the queried k-mers before insertion.
   */
  T contains_insert_thresh(const std::string& seq, const T threshold)
  {
    return contains_insert_thresh(seq.c_str(), seq.size(), threshold);
  }

  /**
   * Get the count of a k-mer and then increment the count if it's not
   * above the threshold.
   *
   * @param hashes Integer array of the k-mer's hash values. Array size
   * should equal the hash_num argument used when the Bloom filter was
   * constructed.
   * @param threshold The threshold.
   *
   * @return The count of the queried k-mer before insertion.
   */
  T contains_insert_thresh(const uint64_t* hashes, const T threshold)
  {
    return counting_bloom_filter.contains_insert_thresh(hashes, threshold);
  }

  /**
   * Get the count of a k-mer and then increment the count if it's not
   * above the threshold.
   *
   * @param hashes Integer vector of the k-mer's hash values.
   * @param threshold The threshold.
   *
   * @return The count of the queried k-mer before insertion.
   */
  T contains_insert_thresh(const std::vector<uint64_t>& hashes,
                           const T threshold)
  {
    return counting_bloom_filter.contains_insert_thresh(hashes, threshold);
  }

  /** Get filter size in bytes. */
  size_t get_bytes() const { return counting_bloom_filter.get_bytes(); }
  /** Get population count, i.e. the number of counters >0 in the filter. */
  uint64_t get_pop_cnt(T threshold = 1) const
  {
    return counting_bloom_filter.get_pop_cnt(threshold);
  }
  /** Get the fraction of the filter occupied by >0 counters. */
  double get_occupancy(T threshold = 1) const
  {
    return counting_bloom_filter.get_occupancy(threshold);
  }
  /** Get the number of hash values per element. */
  unsigned get_hash_num() const { return counting_bloom_filter.get_hash_num(); }
  /** Get the query false positive rate for elements with count >= threshold.
   *
   * @param threshold The threshold.
   *
   * @return The false positive rate.
   */
  double get_fpr(T threshold = 1) const
  {
    return counting_bloom_filter.get_fpr(threshold);
  }
  /** Get the k-mer size used. */
  unsigned get_k() const { return k; }
  /** Get the name of the hash function used. */
  const std::string& get_hash_fn() const
  {
    return counting_bloom_filter.get_hash_fn();
  }
  /** Get a reference to the underlying vanilla Counting Bloom filter. */
  CountingBloomFilter<T>& get_counting_bloom_filter()
  {
    return counting_bloom_filter;
  }

  /**
   * Save the Bloom filter to a file that can be loaded in the future.
   *
   * @param path Filepath to store filter at.
   */
  void save(const std::string& path);

  /**
   * Check whether the file at the given path is a saved k-mer counting Bloom
   * filter.
   *
   * @param path Filepath to check.
   */
  static bool is_bloom_file(const std::string& path)
  {
    return btllib::BloomFilter::check_file_signature(
      path, KMER_COUNTING_BLOOM_FILTER_SIGNATURE);
  }

private:
  KmerCountingBloomFilter(const std::shared_ptr<BloomFilterInitializer>& bfi);

  unsigned k = 0;
  CountingBloomFilter<T> counting_bloom_filter;
};

} // namespace btllib

#include "counting_bloom_filter-inl.hpp"

#endif

#ifndef BTLLIB_MI_BLOOM_FILTER_HPP
#define BTLLIB_MI_BLOOM_FILTER_HPP

#include "nthash.hpp"
#include "status.hpp"

#include "cpptoml.h"

#include <climits>
#include <cstdlib>

#include <sdsl/bit_vector_il.hpp>
#include <sdsl/rank_support.hpp>

namespace btllib {

static const char* const MI_BLOOM_FILTER_SIGNATURE = "[BTLMIBloomFilter_v2]";

static const unsigned PLACEHOLDER_NEWLINES_MIBF = 50;

class MIBloomFilterInitializer
{

  /// @cond HIDDEN_SYMBOLS
public:
  MIBloomFilterInitializer(const std::string& path,
                           const std::string& signature)
    : path(path)
    , ifs_id_arr(path)
    , table(parse_header(signature))
  {
  }

  static bool check_file_signature(std::ifstream& ifs,
                                   const std::string& expected_signature,
                                   std::string& file_signature);

  std::string path;
  std::ifstream ifs_id_arr;
  std::shared_ptr<cpptoml::table> table;

  MIBloomFilterInitializer(const MIBloomFilterInitializer&) = delete;
  MIBloomFilterInitializer(MIBloomFilterInitializer&&) = default;

  MIBloomFilterInitializer& operator=(const MIBloomFilterInitializer&) = delete;
  MIBloomFilterInitializer& operator=(MIBloomFilterInitializer&&) = default;

private:
  /** Parse a Bloom filter file header. Useful for implementing Bloom filter
   * variants. */
  std::shared_ptr<cpptoml::table> parse_header(const std::string& signature);
};
/// @endcond

template<typename T>
class MIBloomFilter
{
public:
  static const T MASK = T(1) << (sizeof(T) * 8 - 1);
  static const T ANTI_MASK = (T)~MASK;

  static const T STRAND = T(1) << (sizeof(T) * 8 - 2);
  static const T ANTI_STRAND = (T)~STRAND;

  static const T ID_MASK = ANTI_STRAND & ANTI_MASK;

  static const unsigned BLOCKSIZE = 512;

  /** Construct a dummy multi-indexed Bloom filter (e.g. as a default argument).
   */
  MIBloomFilter() {}

  /**
   * Construct an empty multi-indexed Bloom filter of given size.
   *
   * @param bv_size Filter size in bytes.
   * @param hash_num Number of hash functions to be used.
   * @param hash_fn Name of the hash function used. Used for metadata. Optional.
   */
  MIBloomFilter(size_t bv_size, unsigned hash_num, std::string hash_fn = "");

  /**
   * Construct a multi-indexed Bloom filter with a prebuilt interleaved bit
   * vector.
   *
   * @param hash_num Number of hash functions to be used.
   * @param hash_fn Name of the hash function used. Used for metadata. Optional.
   */
  MIBloomFilter(sdsl::bit_vector& bit_vector,
                unsigned hash_num,
                std::string hash_fn = "");

  /**
   * Load a multi-indexed Bloom filter from a file.
   *
   * @param path Filepath to load from.
   */
  explicit MIBloomFilter(const std::string& path);

  /**
   * Transform bit vector to interleaved bit vector and create ID array of size
   * of interleaved bit vector. It can be called only once and bit vector
   * insertions cannot be made after calling this function.
   */
  void complete_bv_insertion();

  /**
   * Required to pass to saturation stage.
   */
  void complete_id_insertion() { id_insertion_completed = true; }

  /**
   * Insert an element's hash values to bit vector.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the multi-indexed Bloom filter was constructed.
   */
  void insert_bv(const uint64_t* hashes);

  /**
   * Insert an element's hash values to bit vector.
   *
   * @param hashes Integer vector of hash values. Array size should equal the
   * hash_num argument used when the multi-indexed Bloom filter was constructed.
   */
  void insert_bv(const std::vector<uint64_t>& hashes)
  {
    insert_bv(hashes.data());
  }

  /**
   * Check for presence of an element's hash values in the bit vector
   * This function can run after bit vector insertion is completed.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the multi-indexed Bloom filter was constructed.
   */
  bool bv_contains(const uint64_t* hashes);

  /**
   * Check for presence of an element's hash values in the bit vector
   * This function can run after bit vector insertion is completed.
   *
   * @param hashes Integer vector of hash values. Array size should equal the
   * hash_num argument used when the multi-indexed Bloom filter was constructed.
   */
  bool bv_contains(const std::vector<uint64_t>& hashes)
  {
    return bv_contains(hashes.data());
  }

  /**
   * Insert ID to ID array corresponding to the buckets calculated by hashes.
   * If other ID's were inserted to the same ID bucket previously,
   * even probability for all ID's is assured by random sampling.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   * @param ID The ID that will be inserted into the ID array.
   */
  void insert_id(const uint64_t* hashes, const T& id);

  /**
   * Insert ID to ID array corresponding to the buckets calculated by hashes.
   * If other ID's were inserted to the same ID bucket previously,
   * even probability for all ID's is assured by random sampling.
   *
   * @param hashes Integer vector of hash values.
   * @param ID The ID that will be inserted into the ID array.
   */
  void insert_id(const std::vector<uint64_t>& hashes, const T& id)
  {
    insert_id(hashes.data(), id);
  }

  /**
   * Get the ID's for corresponding to the hashes.
   *
   * @param hashes Integer vector of hash values.
   */
  std::vector<T> get_id(const uint64_t* hashes);

  /**
   * Get the ID's for corresponding to the hashes.
   *
   * @param hashes Integer vector of hash values.
   */
  std::vector<T> get_id(const std::vector<uint64_t>& hashes)
  {
    return get_id(hashes.data());
  }

  /**
   * Inserts saturation if ID is not represented after trying to survive.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * @param ID is the ID to look for.
   */
  void insert_saturation(const uint64_t* hashes, const T& id);

  /**
   * Inserts saturation if ID is not represented after trying to survive.
   *
   * @param hashes Integer vector of hash values.
   * @param ID is the ID to look for.
   */
  void insert_saturation(const std::vector<uint64_t>& hashes, const T& id)
  {
    insert_saturation(hashes.data(), id);
  }

  /**
   * Save the multi-indexed Bloom filter to a file that can be loaded in the
   * future.
   *
   * @param path Filepath to store filter at.
   */
  void save(const std::string& path);

  /** Get population count, i.e. the number of 1 bits in the filter. */
  uint64_t get_pop_cnt();

  /** Get saturation count, i.e. the number of ID's having 1 on leftest-bit in
   * the filter. */
  uint64_t get_pop_saturated_cnt();

  /** Get the number of hash values per element. */
  unsigned get_hash_num() const { return hash_num; }

  /** Get the k-mer size used. */
  unsigned get_k() const { return kmer_size; }

  /** Get the name of the hash function used. */
  const std::string& get_hash_fn() const { return hash_fn; }

  /** Returns the occurence count for each ID in the miBF */
  std::vector<size_t> get_id_occurence_count(const bool& include_saturated);

  /** Returns an a filter size large enough to maintain an occupancy specified
   */
  static size_t calc_optimal_size(size_t entries,
                                  unsigned hash_num,
                                  double occupancy);

private:
  MIBloomFilter(const std::shared_ptr<MIBloomFilterInitializer>& mibfi);
  static void save(const std::string& path,
                   const cpptoml::table& table,
                   const char* data,
                   size_t n);
  std::vector<uint64_t> get_rank_pos(const uint64_t* hashes) const;
  uint64_t get_rank_pos(const uint64_t hash) const
  {
    return bv_rank_support(hash % il_bit_vector.size());
  }
  std::vector<T> get_data(const std::vector<uint64_t>& rank_pos) const;
  T get_data(const uint64_t& rank) const { return id_array[rank]; }
  void set_data(const uint64_t& pos, const T& id);
  void set_saturated(const uint64_t* hashes);

  size_t id_array_size = 0;
  size_t bv_size = 0;
  unsigned kmer_size = 0;
  unsigned hash_num = 0;
  std::string hash_fn;

  sdsl::bit_vector bit_vector;
  sdsl::bit_vector_il<BLOCKSIZE> il_bit_vector;
  sdsl::rank_support_il<1> bv_rank_support;
  std::unique_ptr<std::atomic<uint16_t>[]> counts_array;
  std::unique_ptr<std::atomic<T>[]> id_array;

  bool bv_insertion_completed = false, id_insertion_completed = false;
};

} // namespace btllib

#include "mi_bloom_filter-inl.hpp"

#endif

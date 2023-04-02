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

bool
MIBloomFilterInitializer::check_file_signature(
  std::ifstream& ifs,
  const std::string& expected_signature,
  std::string& file_signature)
{
  std::getline(ifs, file_signature);
  return file_signature == expected_signature;
}

std::shared_ptr<cpptoml::table>
MIBloomFilterInitializer::parse_header(const std::string& expected_signature)
{
  check_file_accessibility(path);
  btllib::check_error(ifs_id_arr.fail(),
                      "MIBloomFilterInitializer: failed to open " + path);

  std::string file_signature;
  if (!check_file_signature(ifs_id_arr, expected_signature, file_signature)) {
    log_error(std::string("File signature does not match (possibly version "
                          "mismatch) for file:\n") +
              path + '\n' + "Expected signature:\t" + expected_signature +
              '\n' + "File signature:    \t" + file_signature);
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }

  /* Read bloom filter line by line until it sees "[HeaderEnd]"
   *   which is used to mark the end of the header section and
   *     assigns the header to a char array*/
  std::string toml_buffer(file_signature + '\n');
  std::string line;
  bool header_end_found = false;
  while (bool(std::getline(ifs_id_arr, line))) {
    toml_buffer.append(line + '\n');
    if (line == "[HeaderEnd]") {
      header_end_found = true;
      break;
    }
  }
  if (!header_end_found) {
    log_error("Pre-built multi-index Bloom filter does not have the correct "
              "header end.");
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
  for (unsigned i = 0; i < PLACEHOLDER_NEWLINES_MIBF; i++) {
    std::getline(ifs_id_arr, line);
  }

  // Send the char array to a stringstream for the cpptoml parser to parse
  std::istringstream toml_stream(toml_buffer);
  cpptoml::parser toml_parser(toml_stream);
  const auto header_config = toml_parser.parse();

  // Obtain header values from toml parser and assign them to class members
  const auto header_string =
    file_signature.substr(1, file_signature.size() - 2); // Remove [ ]
  return header_config->get_table(header_string);
}

template<typename T>
MIBloomFilter<T>::MIBloomFilter(const std::string& path)
  : MIBloomFilter<T>::MIBloomFilter(
      std::make_shared<MIBloomFilterInitializer>(path,
                                                 MI_BLOOM_FILTER_SIGNATURE))
{
}

template<typename T>
inline MIBloomFilter<T>::MIBloomFilter(
  const std::shared_ptr<MIBloomFilterInitializer>& mibfi)
  : id_array_size(
      *(mibfi->table->get_as<decltype(id_array_size)>("id_array_size")))
  , kmer_size(*(mibfi->table->get_as<decltype(kmer_size)>("kmer_size")))
  , hash_num(*(mibfi->table->get_as<decltype(hash_num)>("hash_num")))
  , hash_fn(mibfi->table->contains("hash_fn")
              ? *(mibfi->table->get_as<decltype(hash_fn)>("hash_fn"))
              : "")

  , id_array(new std::atomic<T>[id_array_size])
  , bv_insertion_completed(
      static_cast<bool>(*(mibfi->table->get_as<int>("bv_insertion_completed"))))
  , id_insertion_completed(
      static_cast<bool>(*(mibfi->table->get_as<int>("id_insertion_completed"))))
{
  // read id array
  mibfi->ifs_id_arr.read((char*)id_array.get(),
                         std::streamsize(id_array_size * sizeof(T)));
  // read bv and bv rank support
  sdsl::load_from_file(il_bit_vector, mibfi->path + ".sdsl");
  bv_rank_support = sdsl::rank_support_il<1>(&il_bit_vector);

  // init counts array
  counts_array = std::unique_ptr<std::atomic<uint16_t>[]>(
    new std::atomic<uint16_t>[id_array_size]);
  std::memset(
    (void*)counts_array.get(), 0, id_array_size * sizeof(counts_array[0]));

  log_info(
    "MIBloomFilter: Bit vector size: " + std::to_string(il_bit_vector.size()) +
    "\nPopcount: " + std::to_string(get_pop_cnt()));
}

template<typename T>
inline MIBloomFilter<T>::MIBloomFilter(size_t bv_size,
                                       unsigned hash_num,
                                       std::string hash_fn)
  : bv_size(bv_size)
  , hash_num(hash_num)
  , hash_fn(std::move(hash_fn))
{
  bit_vector = sdsl::bit_vector(bv_size);
}

template<typename T>
inline MIBloomFilter<T>::MIBloomFilter(sdsl::bit_vector& bit_vector,
                                       unsigned hash_num,
                                       std::string hash_fn)
  : bit_vector(bit_vector)
  , hash_num(hash_num)
  , hash_fn(std::move(hash_fn))
{
  complete_bv_insertion();
}

template<typename T>
inline void
MIBloomFilter<T>::insert_bv(const uint64_t* hashes)
{
  assert(!bv_insertion_completed);
  // check array size = hash_num
  for (unsigned i = 0; i < hash_num; ++i) {
    uint64_t pos = hashes[i] % bit_vector.size();
    uint64_t* data_index = bit_vector.data() + (pos >> 6); // NOLINT
    uint64_t bit_mask_value = (uint64_t)1 << (pos & 0x3F); // NOLINT
    (void)(__sync_fetch_and_or(data_index, bit_mask_value) >>
             (pos & 0x3F) & // NOLINT
           1);              // NOLINT
  }
}
template<typename T>
inline bool
MIBloomFilter<T>::bv_contains(const uint64_t* hashes)
{
  assert(bv_insertion_completed);
  for (unsigned i = 0; i < hash_num; i++) {
    uint64_t pos = hashes[i] % il_bit_vector.size();
    if (il_bit_vector[pos] == 0) {
      return false;
    }
  }
  return true;
}
template<typename T>
inline void
MIBloomFilter<T>::complete_bv_insertion()
{
  assert(!id_insertion_completed);
  bv_insertion_completed = true;

  il_bit_vector = sdsl::bit_vector_il<BLOCKSIZE>(bit_vector);
  bv_rank_support = sdsl::rank_support_il<1>(&il_bit_vector);
  id_array_size = get_pop_cnt();
  id_array =
    std::unique_ptr<std::atomic<T>[]>(new std::atomic<T>[id_array_size]);
  std::memset((void*)id_array.get(), 0, id_array_size * sizeof(std::atomic<T>));
  counts_array = std::unique_ptr<std::atomic<uint16_t>[]>(
    new std::atomic<uint16_t>[id_array_size]);
  std::memset(
    (void*)counts_array.get(), 0, id_array_size * sizeof(counts_array[0]));
}
template<typename T>
inline void
MIBloomFilter<T>::insert_id(const uint64_t* hashes, const T& id)
{
  assert(bv_insertion_completed && !id_insertion_completed);

  uint rand = std::rand(); // NOLINT
  for (unsigned i = 0; i < hash_num; ++i) {
    uint64_t rank = get_rank_pos(hashes[i]);
    uint16_t count = ++counts_array[rank];
    T random_num = (rand ^ hashes[i]) % count;
    if (random_num == count - 1) {
      set_data(rank, id);
    }
  }
}
template<typename T>
inline std::vector<T>
MIBloomFilter<T>::get_id(const uint64_t* hashes)
{
  return get_data(get_rank_pos(hashes));
}
template<typename T>
inline void
MIBloomFilter<T>::insert_saturation(const uint64_t* hashes, const T& id)
{
  assert(id_insertion_completed);
  std::vector<uint64_t> rank_pos = get_rank_pos(hashes);
  std::vector<T> results = get_data(rank_pos);
  std::vector<T> replacement_ids(hash_num);
  bool value_found = false;
  std::vector<T> seen_set(hash_num);

  for (unsigned i = 0; i < hash_num; i++) {
    T current_result = results[i] & (btllib::MIBloomFilter<T>::ANTI_MASK &
                                     btllib::MIBloomFilter<T>::ANTI_STRAND);
    // break if ID exists
    if (current_result == id) {
      value_found = true;
      break;
    }
    // if haven't seen before add to seen set
    if (find(seen_set.begin(), seen_set.end(), current_result) ==
        seen_set.end()) {
      seen_set.push_back(current_result);
    }
    // if have seen before add to replacement IDs
    else {
      replacement_ids.push_back(current_result);
    }
  }
  // if value not found try to survive
  if (!value_found) {
    uint64_t replacement_pos = id_array_size;
    T min_count = std::numeric_limits<T>::min();
    for (unsigned i = 0; i < hash_num; i++) {
      T current_result = results[i] & btllib::MIBloomFilter<T>::ANTI_MASK;
      if (find(replacement_ids.begin(),
               replacement_ids.end(),
               current_result) != replacement_ids.end()) {
        if (min_count < counts_array[rank_pos[i]]) {
          min_count = counts_array[rank_pos[i]];
          replacement_pos = rank_pos[i];
        }
      }
    }
    if (replacement_pos != id_array_size) {
      set_data(replacement_pos, id);
      ++counts_array[replacement_pos];
    } else {
      set_saturated(hashes);
    }
  }
}
template<typename T>
inline void
MIBloomFilter<T>::set_data(const uint64_t& pos, const T& id)
{
  T old_value;
  do {
    old_value = id_array[pos];
  } while (!(id_array[pos].compare_exchange_strong(
    old_value, old_value > MASK ? (id | MASK) : id)));
}
template<typename T>
inline void
MIBloomFilter<T>::set_saturated(const uint64_t* hashes)
{
  for (unsigned i = 0; i < hash_num; ++i) {
    uint64_t pos = bv_rank_support(hashes[i] % il_bit_vector.size());
    id_array[pos].fetch_or(MASK);
  }
}
template<typename T>
inline std::vector<uint64_t>
MIBloomFilter<T>::get_rank_pos(const uint64_t* hashes) const
{
  std::vector<uint64_t> rank_pos(hash_num);
  for (unsigned i = 0; i < hash_num; ++i) {
    uint64_t pos = hashes[i] % il_bit_vector.size();
    rank_pos[i] = bv_rank_support(pos);
  }
  return rank_pos;
}
template<typename T>
inline std::vector<T>
MIBloomFilter<T>::get_data(const std::vector<uint64_t>& rank_pos) const
{
  std::vector<T> results(hash_num);
  for (unsigned i = 0; i < hash_num; ++i) {
    results[i] = id_array[rank_pos[i]];
  }
  return results;
}
template<typename T>
inline void
MIBloomFilter<T>::save(const std::string& path,
                       const cpptoml::table& table,
                       const char* data,
                       size_t n)
{
  std::ofstream ofs(path.c_str(), std::ios::out | std::ios::binary);

  ofs << table << "[HeaderEnd]\n";
  for (unsigned i = 0; i < PLACEHOLDER_NEWLINES_MIBF; i++) {
    if (i == 1) {
      ofs << "  <binary data>";
    }
    ofs << '\n';
  }

  ofs.write(data, std::streamsize(n));
}

template<typename T>
inline void
MIBloomFilter<T>::save(const std::string& path)
{
  /* Initialize cpptoml root table
   *     Note: Tables and fields are unordered
   *         Ordering of table is maintained by directing the table
   *             to the output stream immediately after completion  */
  auto root = cpptoml::make_table();

  /* Initialize bloom filter section and insert fields
   *       and output to ostream */
  auto header = cpptoml::make_table();
  header->insert("id_array_size", id_array_size);
  header->insert("hash_num", get_hash_num());
  header->insert("kmer_size", get_k());
  header->insert("bv_insertion_completed",
                 static_cast<int>(bv_insertion_completed));
  header->insert("id_insertion_completed",
                 static_cast<int>(id_insertion_completed));

  if (!hash_fn.empty()) {
    header->insert("hash_fn", get_hash_fn());
  }
  std::string header_string = MI_BLOOM_FILTER_SIGNATURE;
  header_string =
    header_string.substr(1, header_string.size() - 2); // Remove [ ]
  root->insert(header_string, header);
  save(path, *root, (char*)id_array.get(), id_array_size * sizeof(id_array[0]));
  sdsl::store_to_file(il_bit_vector, path + ".sdsl");
}

template<typename T>
inline uint64_t
MIBloomFilter<T>::get_pop_cnt()
{
  assert(bv_insertion_completed);
  size_t index = il_bit_vector.size() - 1;
  while (il_bit_vector[index] == 0) {
    --index;
  }
  return bv_rank_support(index) + 1;
}
template<typename T>
inline uint64_t
MIBloomFilter<T>::get_pop_saturated_cnt()
{
  size_t count = 0;
  for (size_t i = 0; i < id_array_size; ++i) {
    if (id_array[i] >= MASK) {
      ++count;
    }
  }
  return count;
}
template<typename T>
inline std::vector<size_t>
MIBloomFilter<T>::get_id_occurence_count(const bool& include_saturated)
{
  assert(bv_insertion_completed);
  std::vector<size_t> ret_vec(MASK - 1, 0);
#pragma omp parallel for default(none)                                         \
  shared(id_array_size, include_saturated, id_array, ret_vec)
  for (size_t k = 0; k < id_array_size; k++) {
    if (!include_saturated && id_array[k] > ANTI_MASK) {
      continue;
    }
    ret_vec[id_array[k] & ANTI_MASK]++;
  }
  return ret_vec;
}
template<typename T>
inline size_t
MIBloomFilter<T>::calc_optimal_size(size_t entries,
                                    unsigned hash_num,
                                    double occupancy)
{
  auto non_64_approx_val =
    size_t(-double(entries) * double(hash_num) / log(1.0 - occupancy));
  const int magic = 64;
  return non_64_approx_val + (magic - non_64_approx_val % magic);
}
} // namespace btllib
#endif

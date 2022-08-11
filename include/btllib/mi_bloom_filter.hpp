#ifndef BTLLIB_MI_BLOOM_FILTER_HPP
#define BTLLIB_MI_BLOOM_FILTER_HPP

#include "nthash.hpp"
#include "status.hpp"
#include <stdlib.h>

#include <sdsl/bit_vector_il.hpp>
#include <sdsl/rank_support.hpp>

namespace btllib {

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

  /// @cond HIDDEN_SYMBOLS
#pragma pack(1) // to maintain consistent values across platforms
  struct FileHeader
  {
    char magic[8]; // NOLINT
    uint32_t hlen; // header length (including spaced seeds)
    uint64_t size;
    uint32_t nhash;
    uint32_t kmer;
    uint32_t version;
    uint32_t bucket_size;
    uint32_t ID_size;
    uint8_t allowed_miss;
  };
  /// @endcond

  /** Construct a dummy multi-indexed Bloom filter (e.g. as a default argument).
   */
  MIBloomFilter() {}

  /**
   * Construct an empty Bloom filter of given size.
   *
   * @param bv_size Filter size in bv_size.
   * @param hash_fn Name of the hash function used. Used for metadata. Optional.
   */
  MIBloomFilter(size_t bv_size, unsigned hash_num, std::string hash_fn = "");

  /**
  * Construct a Bloom filter with a prebuilt interleaved bit vector.
  *
  * @param bv_size Filter size in bv_size.
  * @param hash_fn Name of the hash function used. Used for metadata. Optional.
  */
  MIBloomFilter(sdsl::bit_vector& bit_vector, unsigned hash_num, std::string hash_fn = "");

  /**
  * Load a Bloom filter from a file.
  *
  * @param path Filepath to load from.
  */
  explicit MIBloomFilter(const std::string& path);

  /**
   * Transform bit vector to interleaved bit vector and create ID array of size
   * of interleaved bit vector. It can be called only once and bit vector
   * insertions cannot be made after calling this function.
   */
  bool complete_bv_insertion();

  /**
   * Insert an element's hash values to bit vector.
   *
   * @param hashes Integer array of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   */
  void insert_bv(const uint64_t* hashes);

  /**
   * Insert an element's hash values to bit vector.
   *
   * @param hashes Integer vector of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
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
   * hash_num argument used when the Bloom filter was constructed.
   */
  bool bv_contains(const uint64_t* hashes);

  /**
   * Check for presence of an element's hash values in the bit vector
   * This function can run after bit vector insertion is completed.
   *
   * @param hashes Integer vector of hash values. Array size should equal the
   * hash_num argument used when the Bloom filter was constructed.
   */
  bool bv_contains(const std::vector<uint64_t>& hashes)
  {
    return bv_contains(hashes.data());
  }

  /**
  * Insert an element's hash values.
  *
  * @param hashes Integer array of hash values. Array size should equal the
  * hash_num argument used when the Bloom filter was constructed.
  */
  void insert_id(const uint64_t* hashes, T& ID);

  /**
  * Insert an element's hash values.
  *
  * @param hashes Integer vector of hash values.
  */
  void insert_id(const std::vector<uint64_t>& hashes, T& ID) { insert_id(hashes.data(), ID); }

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
  std::vector<T> get_id(const std::vector<uint64_t>& hashes) { return get_id(hashes.data());}
 
  /**
  * Save the Bloom filter to a file that can be loaded in the future.
  *
  * @param path Filepath to store filter at.
  */
  void save(const std::string& path);
 
  /** Get population count, i.e. the number of 1 bits in the filter. */
  uint64_t get_pop_cnt();

private:
  std::vector<uint64_t> get_rank_pos(const uint64_t* hashes) const;
  uint64_t get_rank_pos(const uint64_t hash) const { return bv_rank_support(hash % il_bit_vector.size()); }
  std::vector<T> get_data(const std::vector<uint64_t>& rank_pos) const;
  T get_data(uint64_t rank) const { return id_array[rank]; }
  void set_data(uint64_t pos, T id);

  void write_header(std::ofstream& out) const;

  size_t bv_size = 0;
  unsigned kmer_size = 0;
  unsigned hash_num;
  std::string hash_fn;

  sdsl::bit_vector bit_vector;
  sdsl::bit_vector_il<BLOCKSIZE> il_bit_vector;
  sdsl::rank_support_il<1> bv_rank_support;
  std::vector<uint32_t> counts_array;
  T* id_array;

  bool bv_insertion_completed = false, id_insertion_completed = false;
  static const uint32_t MI_BLOOM_FILTER_VERSION = 1;
};

template<typename T>
inline MIBloomFilter<T>::MIBloomFilter(size_t bv_size,
                                       unsigned hash_num,
                                       std::string hash_fn)
  : bv_size(bv_size)
  , hash_num(hash_num)
  , hash_fn(hash_fn)
{
  bit_vector = sdsl::bit_vector(bv_size);
}

template<typename T>
inline MIBloomFilter<T>::MIBloomFilter(sdsl::bit_vector& bit_vector,
							unsigned hash_num,
							std::string hash_fn)
  : bit_vector(bit_vector)
  , hash_num(hash_num)
  , hash_fn(hash_fn)
{
  complete_bv_insertion();
}

template<typename T>
MIBloomFilter<T>::MIBloomFilter(const std::string& filter_file_path)
     // TODO: make more streamlined
  {
#pragma omp parallel for default(none) shared(filter_file_path)
    for (unsigned i = 0; i < 2; ++i) {
      if (i == 0) {
        FILE* file = fopen(filter_file_path.c_str(), "rbe");
        check_error(file == nullptr,
                    "MIBloomFilter: File " + filter_file_path +
                      " could not be read.");

        FileHeader header;
        check_error(fread(&header, sizeof(struct FileHeader), 1, file) != 1,
                    "MIBloomFilter: Failed to load header.");
        log_info("/: Loading header...");

        const int magic_nine = 9;
        char magic[magic_nine];
        const int magic_eight = 8;
        memcpy(magic, header.magic, magic_eight);
        magic[magic_eight] = '\0';

        log_info("MIBloomFilter: Loaded header\nmagic: " + std::string(magic) +
                 "\nhlen: " + std::to_string(header.hlen) +
                 "\nsize: " + std::to_string(header.size) +
                 "\nnhash: " + std::to_string(header.nhash) +
                 "\nkmer: " + std::to_string(header.kmer));
        hash_num = header.nhash;
        kmer_size = header.kmer;
        bv_size = header.size;
        id_array = new T[bv_size]();

	// TOD: Doesnt read spaced seeds!!!
	/*
	check_error(
          header.hlen != (sizeof(FileHeader) + kmer_size * m_sseeds.size()),
          "MIBloomFilter: header length: " + std::to_string(header.hlen) +
            " does not match expected length: " +
            std::to_string(sizeof(FileHeader) + kmer_size * m_sseeds.size()) +
            " (likely version mismatch).");
	*/

        check_error(
          header.hlen != sizeof(FileHeader),
          "MIBloomFilter: header length: " + std::to_string(header.hlen) +
            " does not match expected length: " +
            std::to_string(sizeof(FileHeader)) +
            " (likely version mismatch).");

        check_error(strcmp(magic, "MIBLOOMF") != 0,
                    "MIBloomFilter: Bloom filter type does not matc.");

        check_error(header.version != MI_BLOOM_FILTER_VERSION,
                    "MIBloomFilter: Bloom filter version does not match: " +
                      std::to_string(header.version) + " expected " +
                      std::to_string(MI_BLOOM_FILTER_VERSION) + ".");

        log_info("MIBloomFilter: Loading data vector");

        long int l_cur_pos = ftell(file);
        fseek(file, 0, 2);
        size_t file_size = ftell(file) - header.hlen;
        fseek(file, l_cur_pos, 0);

        check_error(file_size != bv_size * sizeof(T),
                    "MIBloomFilter: " + filter_file_path +
                      " does not match size given by its header. Size: " +
                      std::to_string(file_size) + " vs " +
                      std::to_string(bv_size * sizeof(T)) + " bytes.");

        size_t count_read = fread(id_array, file_size, 1, file);

        check_error(count_read != 1 && fclose(file) != 0,
                    "MIBloomFilter: File " + filter_file_path +
                      " could not be read.");
      }

      else {
        std::string bv_filename = filter_file_path + ".sdsl";
        log_info("MIBloomFilter: Loading sdsl interleaved bit vector from: " +
                 bv_filename);
        load_from_file(il_bit_vector, bv_filename);
        bv_insertion_completed = true;
	bv_rank_support = sdsl::rank_support_il<1>(&il_bit_vector);
  	bv_size = get_pop_cnt();
  	id_array = new T[bv_size]();
  	counts_array = std::vector<uint32_t>(get_pop_cnt(), 0);
      }
    }

    log_info("MIBloomFilter: Bit vector size: " + std::to_string(il_bit_vector.size()) +
             "\nPopcount: " + std::to_string(get_pop_cnt()));
    //m_prob_saturated = pow(double(get_pop_saturated()) / double(get_pop_cnt()), hash_num);
}

template<typename T>
inline void
MIBloomFilter<T>::insert_bv(const uint64_t* hashes)
{
  assert(!bv_insertion_completed);
  // check array size = hash_num
  for (unsigned i = 0; i < hash_num; ++i) {
    uint64_t pos = hashes[i] % bit_vector.size();
    uint64_t* data_index = bit_vector.data() + (pos >> 6);
    uint64_t bit_mask_value = (uint64_t)1 << (pos & 0x3F);
    (void)(__sync_fetch_and_or(data_index, bit_mask_value) >> (pos & 0x3F) & 1);
  }
}
template<typename T>
inline bool
MIBloomFilter<T>::bv_contains(const uint64_t* hashes)
{
  assert(bv_insertion_completed);
  for (unsigned i = 0; i < hash_num; i++) {
    if (il_bit_vector[hashes[i]] == 0) {
      return false;
    }
  }
  return true;
}
template<typename T>
inline bool
MIBloomFilter<T>::complete_bv_insertion()
{
  assert(!id_insertion_completed);
  bv_insertion_completed = true;

  il_bit_vector = sdsl::bit_vector_il<BLOCKSIZE>(bit_vector);
  bv_rank_support = sdsl::rank_support_il<1>(&il_bit_vector);
  bv_size = get_pop_cnt();
  id_array = new T[bv_size]();
  counts_array = std::vector<uint32_t>(get_pop_cnt(), 0);
  return true;
}
template<typename T>
inline void 
MIBloomFilter<T>::insert_id(const uint64_t* hashes, T& ID){
    assert(bv_insertion_completed && !id_insertion_completed);
    //hashSet values;
    //values.set_empty_key(il_bit_vector.size());
    uint rand = std::rand();
    for (unsigned i = 0; i < hash_num; ++i) {
      uint64_t rank = get_rank_pos(hashes[i]);
      uint32_t count = __sync_add_and_fetch(&counts_array[rank], 1);
      T random_num = (rand ^ hashes[i]) % count;
      if (random_num == count - 1) {
        set_data(rank, ID);
      }
    }
}
template<typename T>
inline std::vector<T> MIBloomFilter<T>::get_id(const uint64_t* hashes){
	return get_data(get_rank_pos(hashes));
}
template<typename T>
inline void MIBloomFilter<T>::set_data(uint64_t pos, T ID){
	T old_value;
	do {
		old_value = id_array[pos];
		if (old_value > MASK) { // if the bucket was saturated, insert but keep saturation
			ID |= MASK;
		}
	} while (!__sync_bool_compare_and_swap(&id_array[pos], old_value, ID)); 
}
template<typename T>
inline std::vector<uint64_t> MIBloomFilter<T>::get_rank_pos(const uint64_t* hashes) const{
    std::vector<uint64_t> rank_pos(hash_num);
    for (unsigned i = 0; i < hash_num; ++i) {
      uint64_t pos = hashes[i] % il_bit_vector.size();
      rank_pos[i] = bv_rank_support(pos);
    }
    return rank_pos;
}
template<typename T>
inline std::vector<T> MIBloomFilter<T>::get_data(const std::vector<uint64_t>& rank_pos) const{
    std::vector<T> results(hash_num);
    for (unsigned i = 0; i < hash_num; ++i) {
    	results[i] = id_array[rank_pos[i]];
    }
    return results;
}
template<typename T>
inline void MIBloomFilter<T>::write_header(std::ofstream& out) const
{
    FileHeader header;
    const int magic_num = 8;
    memcpy(header.magic, "MIBLOOMF", magic_num);

    //header.hlen = sizeof(struct FileHeader) + m_kmer_size * m_sseeds.size();
    header.hlen = sizeof(struct FileHeader);
    header.kmer = kmer_size;
    header.size = bv_size;
    header.nhash = hash_num;
    header.version = MI_BLOOM_FILTER_VERSION;


    out.write(reinterpret_cast<char*>(&header), sizeof(struct FileHeader));

    /*for (const auto& s : m_sseeds) {
      out.write(s.c_str(), m_kmer_size);
    }*/
}
/*
 *    * Stores the filter as a binary file to the path specified
 *       * Stores uncompressed because the random data tends to
 *          * compress poorly anyway
 *             */
template<typename T>
inline void MIBloomFilter<T>::save(const std::string& filter_file_path){
//#pragma omp parallel for default(none) shared(filter_file_path)
    for (unsigned i = 0; i < 2; ++i) {
      if (i == 0) {
        std::ofstream my_file(filter_file_path.c_str(),
                              std::ios::out | std::ios::binary);

        assert(my_file);
        write_header(my_file);

        my_file.write(reinterpret_cast<char*>(id_array), bv_size * sizeof(T));

        my_file.close();
        assert(my_file);

        FILE* file = fopen(filter_file_path.c_str(), "rbe");
        check_error(file == nullptr,
                    "MIBloomFilter: " + filter_file_path +
                      " could not be read.");
      } else {
        std::string bv_filename = filter_file_path + ".sdsl";
        store_to_file(il_bit_vector, bv_filename);
      }
    }
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

} // namespace btllib
#endif

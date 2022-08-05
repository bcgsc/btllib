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

  /** Construct a dummy multi-indexed Bloom filter (e.g. as a default argument).
   */
  MIBloomFilter() {}

  /**
   * Construct an empty Bloom filter of given size.
   *
   * @param bytes Filter size in bytes.
   * @param hash_fn Name of the hash function used. Used for metadata. Optional.
   */
  MIBloomFilter(size_t bytes, unsigned hash_num, std::string hash_fn = "");

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
  
  /** Get population count, i.e. the number of 1 bits in the filter. */
  uint64_t get_pop_cnt();

private:
  std::vector<uint64_t> get_rank_pos(const uint64_t* hashes) const;
  uint64_t get_rank_pos(const uint64_t hash) const { return bv_rank_support(hash % il_bit_vector.size()); }
  std::vector<T> get_data(const std::vector<uint64_t>& rank_pos) const;
  T get_data(uint64_t rank) const { return id_array[rank]; }
  void set_data(uint64_t pos, T id);

  size_t bytes = 0;
  unsigned k = 0;
  unsigned hash_num;
  std::string hash_fn;

  sdsl::bit_vector bit_vector;
  sdsl::bit_vector_il<BLOCKSIZE> il_bit_vector;
  sdsl::rank_support_il<1> bv_rank_support;
  std::vector<T> counts_array;
  T* id_array;

  bool bv_insertion_completed = false, id_insertion_completed = false;
};

template<typename T>
inline MIBloomFilter<T>::MIBloomFilter(size_t bytes,
                                       unsigned hash_num,
                                       std::string hash_fn)
  : bytes(bytes)
  , hash_num(hash_num)
  , hash_fn(hash_fn)
{
  bit_vector = sdsl::bit_vector(bytes);
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
  id_array = new T[get_pop_cnt()]();
  counts_array = std::vector<T>(get_pop_cnt(), 0);
  return true;
}
template<typename T>
inline void 
MIBloomFilter<T>::insert_id(const uint64_t* hashes, T& ID){
    assert(bv_insertion_completed && !id_insertion_completed);
    //hashSet values;
    //values.set_empty_key(il_bit_vector.size());
    for (unsigned i = 0; i < hash_num; ++i) {
      uint64_t rank = get_rank_pos(hashes[i]);
      T count = __sync_add_and_fetch(&counts_array[rank], 1);
      T random_num = rand() % count;
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

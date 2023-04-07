#ifndef BTLLIB_AAHASH_HPP
#define BTLLIB_AAHASH_HPP

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace btllib {

class AAHash
{
  static constexpr const char* HASH_FN_NAME = "aahash1";

private:
  /** Initialize internal state of iterator */
  bool init();

  const char* seq;
  size_t seq_len;
  const uint8_t hash_num;
  const uint16_t k;

  size_t pos;
  bool initialized;
  std::unique_ptr<uint64_t[]> hashes_array;

  /**
   * Generate a hash value for the first k-mer in the sequence.
   *
   * @param kmer_seq C array containing the sequence's characters.
   * @param k k-mer size.
   *
   * @return Hash value of k-mer_0.
   */
  uint64_t base_hash(const char* kmer_seq, unsigned k, unsigned level = 1);

  /**
   * Perform a roll operation on the hash value by removing char_out and
   * including char_in.
   *
   * @param hash_value Previous hash value computed for the sequence.
   * @param k k-mer size.
   * @param char_out Character to be removed.
   * @param char_in Character to be included.
   *
   * @return Rolled forward hash value.
   */
  uint64_t roll_forward(uint64_t hash_value,
                        unsigned k,
                        unsigned char char_out,
                        unsigned char char_in,
                        unsigned level = 1);

  /**
   * Extend hash array using a base hash value.
   *
   * @param hash_value Base hash value.
   * @param k k-mer size.
   * @param h Size of the resulting hash array (number of extra hashes minus
   * one).
   * @param hash_array Array of size h for storing the output hashes.
   */
  void extend_hashes(uint64_t hash_value,
                     unsigned k,
                     unsigned h,
                     uint64_t* hash_array);

public:
  /**
   * Constructor.
   * @param seq String of DNA sequence to be hashed.
   * @param hash_num Number of hashes to produce per k-mer.
   * @param k K-mer size.
   * @param pos Position in seq to start hashing from.
   */
  AAHash(const std::string& seq, uint8_t hash_num, uint16_t k, size_t pos = 0)
    : seq(seq.data())
    , seq_len(seq.size())
    , hash_num(hash_num)
    , k(k)
    , pos(pos)
    , initialized(false)
    , hashes_array(new uint64_t[hash_num])
  {
  }

  AAHash(const AAHash& aahash)
    : seq(aahash.seq)
    , seq_len(aahash.seq_len)
    , hash_num(aahash.hash_num)
    , k(aahash.k)
    , initialized(aahash.initialized)
    , hashes_array(new uint64_t[hash_num])
  {
  }

  AAHash(AAHash&&) = default;

  /**
   * Calculate the hash values of current k-mer and advance to the next.
   * AAHash advances one amino acid at a time until it finds a k-mer with valid
   * characters and skips over those with invalid characters.
   * This method must be called before hashes() is accessed, for
   * the first and every subsequent hashed kmer. get_pos() may be called at any
   * time to obtain the position of last hashed k-mer or the k-mer to be hashed
   * if roll() has never been called on this AAHash object. It is important to
   * note that the number of roll() calls is NOT necessarily equal to get_pos(),
   * if there are N or invalid characters in the hashed sequence.
   *
   * @return true on success and false otherwise.
   */
  bool roll();

  const uint64_t* hashes() const { return hashes_array.get(); }
  size_t get_pos() const { return pos; }
  unsigned get_hash_num() const { return hash_num; }
  unsigned get_k() const { return k; }
  uint64_t get_forward_hash() const { return hashes_array[0]; }
};

} // namespace btllib

#endif
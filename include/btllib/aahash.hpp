#ifndef BTLLIB_AAHASH_HPP
#define BTLLIB_AAHASH_HPP

#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "aahash_consts.hpp"
#include "nthash_lowlevel.hpp"

namespace btllib {

inline uint64_t
aahash_base(const char* kmer_seq, unsigned k, unsigned level = 1)
{
  uint64_t hash_value = 0;
  for (unsigned i = 0; i < k; i++) {
    hash_value = srol(hash_value);
    hash_value ^= LEVEL_X_AA_SEED_TABLE[level][(unsigned char)kmer_seq[i]];
  }
  return hash_value;
}

inline uint64_t
aahash_roll(uint64_t hash_value,
            unsigned k,
            unsigned char char_out,
            unsigned char char_in,
            unsigned level = 1)
{
  hash_value = srol(hash_value);
  hash_value ^= LEVEL_X_AA_SEED_TABLE[level][char_in];
  hash_value ^= AA_ROLL_TABLE(char_out, level, k);
  return hash_value;
}

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
  unsigned level;

  size_t pos;
  bool initialized;
  std::unique_ptr<uint64_t[]> hashes_array;

public:
  /**
   * Constructor.
   * @param seq String of DNA sequence to be hashed.
   * @param hash_num Number of hashes to produce per k-mer.
   * @param k K-mer size.
   * @param pos Position in seq to start hashing from.
   * @param level seed level to generate hash.
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
    level = 1;
  }

  AAHash(const std::string& seq,
         uint8_t hash_num,
         uint16_t k,
         unsigned level,
         size_t pos = 0)
    : seq(seq.data())
    , seq_len(seq.size())
    , hash_num(hash_num)
    , k(k)
    , level(level)
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
    std::memcpy(hashes_array.get(),
                aahash.hashes_array.get(),
                hash_num * sizeof(uint64_t));
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
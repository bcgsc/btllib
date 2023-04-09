#include "btllib/aahash.hpp"
#include "btllib/aahash_consts.hpp"
#include "btllib/nthash_lowlevel.hpp"

namespace btllib {

uint64_t
AAHash::base_hash(const char* kmer_seq, unsigned k, unsigned level)
{
  uint64_t hash_value = 0;
  for (unsigned i = 0; i < k; i++) {
    hash_value = srol(hash_value);
    hash_value ^= LEVEL_X_AA_SEED_TABLE[level][(unsigned char)kmer_seq[i]];
  }
  return hash_value;
}

uint64_t
AAHash::roll_forward(uint64_t hash_value,
                     unsigned k,
                     unsigned char char_out,
                     unsigned char char_in,
                     unsigned level)
{
  hash_value = srol(hash_value);
  hash_value ^= LEVEL_X_AA_SEED_TABLE[level][char_in];
  hash_value ^= AA_ROLL_TABLE(char_out, level, k);
  return hash_value;
}

void
AAHash::extend_hashes(uint64_t hash_value,
                      unsigned k,
                      unsigned h,
                      uint64_t* hash_array)
{
  uint64_t temp_value;
  hash_array[0] = hash_value;
  for (unsigned i = 1; i < h; i++) {
    temp_value = hash_value * (i ^ k * MULTISEED);
    temp_value ^= temp_value >> MULTISHIFT;
    hash_array[i] = temp_value;
  }
}

bool
AAHash::init()
{
  if (k > seq_len) {
    pos = std::numeric_limits<std::size_t>::max();
    return false;
  }
  if (pos > seq_len - k) {
    pos = std::numeric_limits<std::size_t>::max();
    return false;
  }
  uint64_t hash_value = base_hash(seq + pos, k, level);
  extend_hashes(hash_value, k, hash_num, hashes_array.get());
  initialized = true;
  return true;
}

bool
AAHash::roll()
{
  if (!initialized) {
    return init();
  }
  if (pos >= seq_len - k) {
    return false;
  }
  if (AA_SEED_Table[(unsigned char)(seq[pos + k])] == AA_SEED__) {
    pos += k;
    return init();
  }
  uint64_t hash_value =
    roll_forward(hashes_array[0], k, seq[pos], seq[pos + k], level);
  extend_hashes(hash_value, k, hash_num, hashes_array.get());
  ++pos;
  return true;
}

}
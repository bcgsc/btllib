#include "btllib/aahash.hpp"
#include "btllib/aahash_consts.hpp"
#include "btllib/nthash_lowlevel.hpp"

namespace btllib {

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
  uint64_t hash_value = aahash_base(seq + pos, k, level);
  nte64(hash_value, k, hash_num, hashes_array.get());
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
  uint64_t hash_value = aahash_roll(hashes_array[0], k, seq[pos], seq[pos + k], level);
  nte64(hash_value, k, hash_num, hashes_array.get());
  ++pos;
  return true;
}

}
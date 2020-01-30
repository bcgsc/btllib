#ifndef BTLLIB_ROLLING_HASH_HPP
#define BTLLIB_ROLLING_HASH_HPP

#include "nthash.hpp"

#include <limits>
#include <string>
#include <vector>

namespace btllib {

class RollingHash;
class SeedRollingHash;
using SpacedSeed = std::vector<unsigned>;
static std::vector<SpacedSeed>
parse_seeds(const std::vector<std::string>& seed_strings);

/**
 * Iterate over hash values for k-mers in a
 * given DNA sequence.
 *
 * This implementation uses ntHash
 * function to efficiently calculate
 * hash values for successive k-mers.
 */
class RollingHash
{

public:
  /**
   * Constructor.
   * @param seq DNA sequence to be hashed
   * @param seq_len length of seq
   * @param k k-mer size
   * @param hash_num number of hashes
   */
  RollingHash(const char* seq, size_t seq_len, unsigned k, unsigned hash_num);

  /**
   * Constructor.
   * @param seq DNA sequence to be hashed
   * @param k k-mer size
   * @param hash_num number of hashes
   */
  RollingHash(const std::string& seq, unsigned k, unsigned hash_num);

  bool roll();

  const uint64_t* hashes() const;

  size_t get_pos() const { return pos; }
  unsigned get_k() const { return k; }
  unsigned get_hash_num() const { return hash_num; }

protected:
  /** Initialize internal state of iterator */
  bool init();

  const char* seq;
  const size_t seq_len;
  const unsigned k;
  const unsigned hash_num;
  size_t pos = 0;
  std::vector<uint64_t> hashes_vector;
  uint64_t forward_hash = 0;
  uint64_t reverse_hash = 0;
};

class SeedRollingHash : public RollingHash
{

public:
  SeedRollingHash(const char* seq,
                  size_t seq_len,
                  unsigned k,
                  const std::vector<SpacedSeed>& seeds,
                  unsigned hash_num_per_seed);
  SeedRollingHash(const std::string& seq,
                  unsigned k,
                  const std::vector<SpacedSeed>& seeds,
                  unsigned hash_num_per_seed);
  SeedRollingHash(const char* seq,
                  size_t seq_len,
                  unsigned k,
                  const std::vector<std::string>& seeds,
                  unsigned hash_num_per_seed);
  SeedRollingHash(const std::string& seq,
                  unsigned k,
                  const std::vector<std::string>& seeds,
                  unsigned hash_num_per_seed);

  unsigned get_hash_num_per_seed() const { return hash_num_per_seed; }

  bool roll();

private:
  bool init();

  const unsigned hash_num_per_seed;
  std::vector<SpacedSeed> seeds;
};

inline RollingHash::RollingHash(const char* seq,
                                size_t seq_len,
                                unsigned k,
                                unsigned hash_num)
  : seq(seq)
  , seq_len(seq_len)
  , k(k)
  , hash_num(hash_num)
{
  hashes_vector.resize(hash_num);
}

inline RollingHash::RollingHash(const std::string& seq,
                                unsigned k,
                                unsigned hash_num)
  : RollingHash(seq.c_str(), seq.size(), k, hash_num)
{}

inline SeedRollingHash::SeedRollingHash(const char* seq,
                                        size_t seq_len,
                                        unsigned k,
                                        const std::vector<SpacedSeed>& seeds,
                                        unsigned hash_num_per_seed)
  : RollingHash(seq, seq_len, k, seeds.size() * hash_num_per_seed)
  , hash_num_per_seed(hash_num_per_seed)
  , seeds(seeds)
{}

inline SeedRollingHash::SeedRollingHash(const std::string& seq,
                                        unsigned k,
                                        const std::vector<SpacedSeed>& seeds,
                                        unsigned hash_num_per_seed)
  : RollingHash(seq, k, seeds.size() * hash_num_per_seed)
  , hash_num_per_seed(hash_num_per_seed)
  , seeds(seeds)
{}

inline SeedRollingHash::SeedRollingHash(const char* seq,
                                        size_t seq_len,
                                        unsigned k,
                                        const std::vector<std::string>& seeds,
                                        unsigned hash_num_per_seed)
  : RollingHash(seq, seq_len, k, seeds.size() * hash_num_per_seed)
  , hash_num_per_seed(hash_num_per_seed)
  , seeds(parse_seeds(seeds))
{}

inline SeedRollingHash::SeedRollingHash(const std::string& seq,
                                        unsigned k,
                                        const std::vector<std::string>& seeds,
                                        unsigned hash_num_per_seed)
  : RollingHash(seq, k, seeds.size() * hash_num_per_seed)
  , hash_num_per_seed(hash_num_per_seed)
  , seeds(parse_seeds(seeds))
{}

static std::vector<SpacedSeed>
parse_seeds(const std::vector<std::string>& seed_strings)
{
  std::vector<SpacedSeed> seed_set;
  for (const auto& seed_string : seed_strings) {
    SpacedSeed seed;
    size_t pos = 0;
    for (const auto& c : seed_string) {
      if (c != '1') {
        seed.push_back(pos);
      }
      ++pos;
    }
    seed_set.push_back(seed);
  }
  return seed_set;
}

// NOLINTNEXTLINE
#define ROLLING_HASH_INIT(CLASS, NTHASH_CALL)                                  \
  inline bool CLASS::init()                                                    \
  {                                                                            \
    if (k > seq_len) {                                                         \
      pos = std::numeric_limits<std::size_t>::max();                           \
      return false;                                                            \
    }                                                                          \
    unsigned posN = 0;                                                         \
    while ((pos < seq_len - k + 1) && !(NTHASH_CALL)) {                        \
      pos += posN + 1;                                                         \
    }                                                                          \
    if (pos > seq_len - k) {                                                   \
      pos = std::numeric_limits<std::size_t>::max();                           \
      return false;                                                            \
    }                                                                          \
    ++pos;                                                                     \
    return true;                                                               \
  }

// NOLINTNEXTLINE
#define ROLLING_HASH_ROLL(CLASS, NTHASH_CALL)                                  \
  inline bool CLASS::roll()                                                    \
  {                                                                            \
    if (pos == 0) {                                                            \
      return init();                                                           \
    }                                                                          \
    if (pos > seq_len - k) {                                                   \
      return false;                                                            \
    }                                                                          \
    if (seed_tab[(unsigned char)(seq[pos + k - 1])] == seedN) {                \
      pos += k;                                                                \
      return init();                                                           \
    }                                                                          \
    (NTHASH_CALL);                                                             \
    ++pos;                                                                     \
    return true;                                                               \
  }

ROLLING_HASH_INIT(RollingHash,
                  NTMC64(seq + pos,
                         k,
                         hash_num,
                         forward_hash,
                         reverse_hash,
                         posN,
                         hashes_vector.data()))
ROLLING_HASH_ROLL(RollingHash,
                  NTMC64(seq[pos - 1],
                         seq[pos - 1 + k],
                         k,
                         hash_num,
                         forward_hash,
                         reverse_hash,
                         hashes_vector.data()))

ROLLING_HASH_INIT(SeedRollingHash,
                  NTMSM64(seq + pos,
                          seeds,
                          k,
                          seeds.size(),
                          hash_num_per_seed,
                          forward_hash,
                          reverse_hash,
                          posN,
                          hashes_vector.data()))
ROLLING_HASH_ROLL(SeedRollingHash,
                  NTMSM64(seq + pos,
                          seeds,
                          seq[pos - 1],
                          seq[pos - 1 + k],
                          k,
                          seeds.size(),
                          hash_num_per_seed,
                          forward_hash,
                          reverse_hash,
                          hashes_vector.data()))

#undef ROLLING_HASH_INIT
#undef ROLLING_HASH_ROLL

inline const uint64_t*
RollingHash::hashes() const
{
  return hashes_vector.data();
}

} // namespace btllib

#endif

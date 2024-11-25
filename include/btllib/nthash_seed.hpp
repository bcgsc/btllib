#pragma once

#include <array>
#include <cstdint>
#include <cstring>
#include <deque>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <btllib/hashing_internals.hpp>
#include <btllib/status.hpp>

namespace btllib::hashing_internals {

using SpacedSeed = std::vector<unsigned>;
using SpacedSeedBlocks = std::vector<std::array<unsigned, 2>>;
using SpacedSeedMonomers = std::vector<unsigned>;

inline void
get_blocks(const std::vector<std::string>& seed_strings,
           std::vector<SpacedSeedBlocks>& blocks,
           std::vector<SpacedSeedMonomers>& monomers)
{
  for (const auto& seed_string : seed_strings) {
    const char pad = seed_string[seed_string.length() - 1] == '1' ? '0' : '1';
    const std::string padded_string = seed_string + pad;
    SpacedSeedBlocks care_blocks, ignore_blocks;
    std::vector<unsigned> care_monos, ignore_monos;
    unsigned i_start = 0;
    bool is_care_block = padded_string[0] == '1';
    for (unsigned pos = 0; pos < padded_string.length(); pos++) {
      if (is_care_block && padded_string[pos] == '0') {
        if (pos - i_start == 1) {
          care_monos.push_back(i_start);
        } else {
          const std::array<unsigned, 2> block{ { i_start, pos } };
          care_blocks.push_back(block);
        }
        i_start = pos;
        is_care_block = false;
      } else if (!is_care_block && padded_string[pos] == '1') {
        if (pos - i_start == 1) {
          ignore_monos.push_back(i_start);
        } else {
          const std::array<unsigned, 2> block{ { i_start, pos } };
          ignore_blocks.push_back(block);
        }
        i_start = pos;
        is_care_block = true;
      }
    }
    const unsigned num_cares = care_blocks.size() * 2 + care_monos.size();
    const unsigned num_ignores =
      ignore_blocks.size() * 2 + ignore_monos.size() + 2;
    if (num_ignores < num_cares) {
      const unsigned string_end = seed_string.length();
      const std::array<unsigned, 2> block{ { 0, string_end } };
      ignore_blocks.push_back(block);
      blocks.push_back(ignore_blocks);
      monomers.push_back(ignore_monos);
    } else {
      blocks.push_back(care_blocks);
      monomers.push_back(care_monos);
    }
  }
}

inline void
parsed_seeds_to_blocks(const std::vector<std::vector<unsigned>>& seeds,
                       unsigned k,
                       std::vector<SpacedSeedBlocks>& blocks,
                       std::vector<SpacedSeedMonomers>& monomers)
{
  std::vector<std::string> seed_strings;
  for (const auto& seed : seeds) {
    std::string seed_string(k, '1');
    for (const auto& i : seed) {
      seed_string[i] = '0';
    }
    seed_strings.push_back(seed_string);
  }
  get_blocks(seed_strings, blocks, monomers);
}

inline void
check_seeds(const std::vector<std::string>& seeds, unsigned k)
{
  for (const auto& seed : seeds) {
    btllib::check_error(seed.length() != k,
                        "SeedNtHash: Spaced seed string length (" +
                          std::to_string(seed.length()) + ") not equal to k=" +
                          std::to_string(k) + " in " + seed);
    const std::string reversed(seed.rbegin(), seed.rend());
    btllib::check_warning(
      seed != reversed,
      "SeedNtHash: Seed " + seed +
        " is not symmetric, reverse-complement hashing will be inconsistent");
  }
}

/**
 * Generate multiple hash values for the input spaced seeds and first k-mer.
 *
 * @param kmer_seq Array of characters representing the k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Container for the forward hash values before including the
 * size-one blocks.
 * @param rh_nomonos Container for the reverse hash values before including the
 * size-one blocks.
 * @param fh_val Container for the forward hash values after including the
 * size-one blocks.
 * @param rh_val Container for the reverse hash values after including the
 * size-one blocks.
 * @param loc_n Location of the first unknown character in the first sequence.
 * @param h_val Array of size m * m2 for storing the output hash values.
 *
 * @return true if all the care positions of the first k-mer are valid,
 * otherwise false.
 */
inline bool
ntmsm64(const char* kmer_seq,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        unsigned k,
        unsigned m,
        unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        unsigned& loc_n,
        uint64_t* h_val)
{
  unsigned i_base;
  uint64_t fh_seed, rh_seed;
  for (unsigned i_seed = 0; i_seed < m; i_seed++) {
    fh_seed = 0;
    rh_seed = 0;
    for (const auto& block : seeds_blocks[i_seed]) {
      uint8_t fh_loc, rh_loc, d;
      uint64_t x;
      unsigned i = block[0];
      unsigned j = block[1];
      switch (j - i) {
        case 2: {
          fh_loc = (CONVERT_TAB[(unsigned char)kmer_seq[i]] << 2) | // NOLINT
                   (CONVERT_TAB[(unsigned char)kmer_seq[i + 1]]);   // NOLINT
          rh_loc =
            (RC_CONVERT_TAB[(unsigned char)kmer_seq[i + 1]] << 2) | // NOLINT
            (RC_CONVERT_TAB[(unsigned char)kmer_seq[i]]);           // NOLINT
          x = DIMER_TAB[fh_loc]; // cppcheck-suppress arrayIndexOutOfBounds
          d = k - i - 2;         // NOLINT
          fh_seed ^= d > 0 ? srol(x, d) : x;
          x = DIMER_TAB[rh_loc]; // cppcheck-suppress arrayIndexOutOfBounds
          d = i;
          rh_seed ^= d > 0 ? srol(x, d) : x;
        } break;
        case 3: {
          fh_loc =
            (CONVERT_TAB[(unsigned char)kmer_seq[i]] << 4) |     // NOLINT
            (CONVERT_TAB[(unsigned char)kmer_seq[i + 1]] << 2) | // NOLINT
            (CONVERT_TAB[(unsigned char)kmer_seq[i + 2]]);       // NOLINT
          rh_loc =
            (RC_CONVERT_TAB[(unsigned char)kmer_seq[i + 2]] << 4) | // NOLINT
            (RC_CONVERT_TAB[(unsigned char)kmer_seq[i + 1]] << 2) | // NOLINT
            (RC_CONVERT_TAB[(unsigned char)kmer_seq[i]]);
          x = TRIMER_TAB[fh_loc]; // cppcheck-suppress arrayIndexOutOfBounds
          d = k - i - 3;          // NOLINT
          fh_seed ^= d > 0 ? srol(x, d) : x;
          x = TRIMER_TAB[rh_loc]; // cppcheck-suppress arrayIndexOutOfBounds
          d = i;
          rh_seed ^= d > 0 ? srol(x, d) : x;
        } break;
        case 4: {
          fh_loc =
            (CONVERT_TAB[(unsigned char)kmer_seq[i]] << 6) |     // NOLINT
            (CONVERT_TAB[(unsigned char)kmer_seq[i + 1]] << 4) | // NOLINT
            (CONVERT_TAB[(unsigned char)kmer_seq[i + 2]] << 2) | // NOLINT
            (CONVERT_TAB[(unsigned char)kmer_seq[i + 3]]);       // NOLINT
          rh_loc =
            (RC_CONVERT_TAB[(unsigned char)kmer_seq[i + 3]] << 6) | // NOLINT
            (RC_CONVERT_TAB[(unsigned char)kmer_seq[i + 2]] << 4) | // NOLINT
            (RC_CONVERT_TAB[(unsigned char)kmer_seq[i + 1]] << 2) | // NOLINT
            (RC_CONVERT_TAB[(unsigned char)kmer_seq[i]]);
          x = TETRAMER_TAB[fh_loc]; // cppcheck-suppress arrayIndexOutOfBounds
          d = k - i - 4;            // NOLINT
          fh_seed ^= d > 0 ? srol(x, d) : x;
          x = TETRAMER_TAB[rh_loc]; // cppcheck-suppress arrayIndexOutOfBounds
          d = i;
          rh_seed ^= d > 0 ? srol(x, d) : x;
        } break;
        default: {
          for (unsigned pos = block[0]; pos < block[1]; pos++) {
            if (kmer_seq[pos] == SEED_N) {
              loc_n = pos;
              return false;
            }
            fh_seed ^= srol_table((unsigned char)kmer_seq[pos], k - 1 - pos);
            rh_seed ^= srol_table((unsigned char)kmer_seq[pos] & CP_OFF, pos);
          }
        }
      }
    }
    fh_nomonos[i_seed] = fh_seed;
    rh_nomonos[i_seed] = rh_seed;
    for (const auto& pos : seeds_monomers[i_seed]) {
      fh_seed ^= srol_table((unsigned char)kmer_seq[pos], k - 1 - pos);
      rh_seed ^= srol_table((unsigned char)kmer_seq[pos] & CP_OFF, pos);
    }
    fh_val[i_seed] = fh_seed;
    rh_val[i_seed] = rh_seed;
    i_base = i_seed * m2;
    h_val[i_base] = canonical(fh_seed, rh_seed);
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {
      h_val[i_base + i_hash] = h_val[i_base] * (i_hash ^ k * MULTISEED);
      h_val[i_base + i_hash] ^= h_val[i_base + i_hash] >> MULTISHIFT;
    }
  }
  return true;
}

#define NTMSM64(ROL_HANDLING, IN_HANDLING, OUT_HANDLING, ROR_HANDLING)         \
  unsigned char char_out, char_in;                                             \
  uint64_t fh_seed, rh_seed;                                                   \
  unsigned i_out, i_in, i_base;                                                \
  for (unsigned i_seed = 0; i_seed < m; i_seed++) {                            \
    ROL_HANDLING /* NOLINT(bugprone-macro-parentheses) */                      \
      for (const auto& block : seeds_blocks[i_seed])                           \
    {                                                                          \
      IN_HANDLING                                                              \
      OUT_HANDLING                                                             \
      fh_seed ^= srol_table(char_out, k - i_out);                              \
      fh_seed ^= srol_table(char_in, k - i_in);                                \
      rh_seed ^= srol_table(char_out & CP_OFF, i_out);                         \
      rh_seed ^= srol_table(char_in & CP_OFF, i_in);                           \
    }                                                                          \
    ROR_HANDLING /* NOLINT(bugprone-macro-parentheses) */                      \
      fh_nomonos[i_seed] = fh_seed;                                            \
    rh_nomonos[i_seed] = rh_seed;                                              \
    for (const auto& pos : seeds_monomers[i_seed]) {                           \
      fh_seed ^= srol_table((unsigned char)kmer_seq[pos + 1], k - 1 - pos);    \
      rh_seed ^= srol_table((unsigned char)kmer_seq[pos + 1] & CP_OFF, pos);   \
    }                                                                          \
    fh_val[i_seed] = fh_seed;                                                  \
    rh_val[i_seed] = rh_seed;                                                  \
    i_base = i_seed * m2;                                                      \
    h_val[i_base] = canonical(fh_seed, rh_seed);                               \
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {                         \
      h_val[i_base + i_hash] = h_val[i_base] * (i_hash ^ k * MULTISEED);       \
      h_val[i_base + i_hash] ^= h_val[i_base + i_hash] >> MULTISHIFT;          \
    }                                                                          \
  }

/**
 * Generate multiple hash values for the input spaced seeds and the next
 * k-mer by performing a forward roll operation.
 *
 * @param kmer_seq Array of characters representing the previous k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Previous forward hash values before including the size-one
 * blocks.
 * @param rh_nomonos Previous reverse hash values before including the size-one
 * blocks.
 * @param fh_val Previous forward hash values after including the size-one
 * blocks.
 * @param rh_val Previous reverse hash values after including the size-one
 * blocks.
 * @param h_val Array of size m * m2 for storing the output hash values.
 */
inline void
ntmsm64(const char* kmer_seq,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        unsigned k,
        unsigned m,
        unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        uint64_t* h_val)
{
  NTMSM64(fh_seed = srol(fh_nomonos[i_seed]); rh_seed = rh_nomonos[i_seed];
          , i_in = block[1];
          char_in = (unsigned char)kmer_seq[i_in];
          , i_out = block[0];
          char_out = (unsigned char)kmer_seq[i_out];
          , rh_seed = sror(rh_seed);)
}

inline void
ntmsm64(const std::deque<char>& kmer_seq,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        unsigned k,
        unsigned m,
        unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        uint64_t* h_val)
{
  NTMSM64(fh_seed = srol(fh_nomonos[i_seed]); rh_seed = rh_nomonos[i_seed];
          , i_in = block[1];
          char_in = (unsigned char)kmer_seq[i_in];
          , i_out = block[0];
          char_out = (unsigned char)kmer_seq[i_out];
          , rh_seed = sror(rh_seed);)
}

/**
 * Generate multiple hash values for the input spaced seeds and the next
 * k-mer by performing a backward roll operation.
 *
 * @param kmer_seq Array of characters representing the previous k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Previous forward hash values before including the size-one
 * blocks.
 * @param rh_nomonos Previous reverse hash values before including the size-one
 * blocks.
 * @param fh_val Previous forward hash values after including the size-one
 * blocks.
 * @param rh_val Previous reverse hash values after including the size-one
 * blocks.
 * @param h_val Array of size m * m2 for storing the output hash values.
 */
inline void
ntmsm64l(const char* kmer_seq,
         const std::vector<SpacedSeedBlocks>& seeds_blocks,
         const std::vector<SpacedSeedMonomers>& seeds_monomers,
         unsigned k,
         unsigned m,
         unsigned m2,
         uint64_t* fh_nomonos,
         uint64_t* rh_nomonos,
         uint64_t* fh_val,
         uint64_t* rh_val,
         uint64_t* h_val)
{
  NTMSM64(fh_seed = fh_nomonos[i_seed]; rh_seed = srol(rh_nomonos[i_seed]);
          , i_in = block[0];
          char_in = (unsigned char)kmer_seq[i_in];
          , i_out = block[1];
          char_out = (unsigned char)kmer_seq[i_out];
          , fh_seed = sror(fh_seed);)
}

inline void
ntmsm64l(const std::deque<char>& kmer_seq,
         const std::vector<SpacedSeedBlocks>& seeds_blocks,
         const std::vector<SpacedSeedMonomers>& seeds_monomers,
         unsigned k,
         unsigned m,
         unsigned m2,
         uint64_t* fh_nomonos,
         uint64_t* rh_nomonos,
         uint64_t* fh_val,
         uint64_t* rh_val,
         uint64_t* h_val)
{
  NTMSM64(fh_seed = fh_nomonos[i_seed]; rh_seed = srol(rh_nomonos[i_seed]);
          , i_in = block[0];
          char_in = (unsigned char)kmer_seq[i_in];
          , i_out = block[1];
          char_out = (unsigned char)kmer_seq[i_out];
          , fh_seed = sror(fh_seed);)
}

/**
 * Generate multiple hash values for the input spaced seeds and the next
 * k-mer by performing a forward peek operation.
 *
 * @param kmer_seq Array of characters representing the previous k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Previous forward hash values before including the size-one
 * blocks.
 * @param rh_nomonos Previous reverse hash values before including the size-one
 * blocks.
 * @param fh_val Previous forward hash values after including the size-one
 * blocks.
 * @param rh_val Previous reverse hash values after including the size-one
 * blocks.
 * @param h_val Array of size m * m2 for storing the output hash values.
 */
inline void
ntmsm64(const char* kmer_seq,
        char in,
        const std::vector<SpacedSeedBlocks>& seeds_blocks,
        const std::vector<SpacedSeedMonomers>& seeds_monomers,
        unsigned k,
        unsigned m,
        unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        uint64_t* h_val)
{
  NTMSM64(
    fh_seed = srol(fh_nomonos[i_seed]); rh_seed = rh_nomonos[i_seed];
    , i_in = block[1];
    if (i_in > k - 1) { char_in = in; } else {
      char_in = (unsigned char)kmer_seq[i_in];
    },
    i_out = block[0];
    char_out = (unsigned char)kmer_seq[i_out];
    , rh_seed = sror(rh_seed);)
}

/**
 * Generate multiple hash values for the input spaced seeds and the next
 * k-mer by performing a backwards peek operation.
 *
 * @param kmer_seq Array of characters representing the previous k-mer.
 * @param seed_seq Array of SpacedSeed objects representing the seeds' blocks.
 * @param monomers List of the positions that represent blocks of size one for
 * each seed.
 * @param k k-mer size.
 * @param m Number of spaced seeds.
 * @param m2 Number of hashes per seed.
 * @param fh_nomonos Previous forward hash values before including the size-one
 * blocks.
 * @param rh_nomonos Previous reverse hash values before including the size-one
 * blocks.
 * @param fh_val Previous forward hash values after including the size-one
 * blocks.
 * @param rh_val Previous reverse hash values after including the size-one
 * blocks.
 * @param h_val Array of size m * m2 for storing the output hash values.
 */
inline void
ntmsm64l(const char* kmer_seq,
         char in,
         const std::vector<SpacedSeedBlocks>& seeds_blocks,
         const std::vector<SpacedSeedMonomers>& seeds_monomers,
         unsigned k,
         unsigned m,
         unsigned m2,
         uint64_t* fh_nomonos,
         uint64_t* rh_nomonos,
         uint64_t* fh_val,
         uint64_t* rh_val,
         uint64_t* h_val)
{
  NTMSM64(
    fh_seed = fh_nomonos[i_seed]; rh_seed = srol(rh_nomonos[i_seed]);
    , i_in = block[0];
    if (i_in > k - 1) { char_in = in; } else {
      char_in = (unsigned char)kmer_seq[i_in];
    },
    i_out = block[1];
    char_out = (unsigned char)kmer_seq[i_out];
    , fh_seed = sror(fh_seed);)
}

} // namespace btllib::hashing_internals

namespace btllib {

using hashing_internals::check_seeds;
using hashing_internals::get_blocks;
using hashing_internals::ntmsm64;
using hashing_internals::ntmsm64l;
using hashing_internals::parsed_seeds_to_blocks;
using hashing_internals::SEED_N;
using hashing_internals::SEED_TAB;

/**
 * Parse each spaced seed pattern into lists of don't care positions. Legacy
 * function used in btllib Bloom filters.
 */
inline std::vector<std::vector<unsigned>>
parse_seeds(const std::vector<std::string>& seed_strings)
{
  std::vector<std::vector<unsigned>> seed_set;
  for (const auto& seed_string : seed_strings) {
    std::vector<unsigned> seed;
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

/**
 * Spaced seed hashing.
 */
class SeedNtHash
{

public:
  /**
   * Construct an ntHash object for spaced seeds.
   * @param seq C-string of the sequence to be hashed
   * @param seq_len Length of the sequence
   * @param seeds Vector of spaced seed patterns as strings (1s as cares, 0s
   * as don't cares, must be of size \p k)
   * @param num_hashes_per_seed Number of hashes to generate per seed
   * @param k K-mer size
   * @param pos Position in seq to start hashing from
   */
  SeedNtHash(const char* seq,
             size_t seq_len,
             const std::vector<std::string>& seeds,
             hashing_internals::NUM_HASHES_TYPE num_hashes_per_seed,
             hashing_internals::K_TYPE k,
             size_t pos = 0)
    : seq(seq)
    , seq_len(seq_len)
    , num_hashes_per_seed(num_hashes_per_seed)
    , k(k)
    , pos(pos)
    , initialized(false)
    , fwd_hash_nomonos(new uint64_t[seeds.size()])
    , rev_hash_nomonos(new uint64_t[seeds.size()])
    , fwd_hash(new uint64_t[seeds.size()])
    , rev_hash(new uint64_t[seeds.size()])
    , hash_arr(new uint64_t[num_hashes_per_seed * seeds.size()])
  {
    check_seeds(seeds, k);
    check_error(seeds[0].size() != k,
                "SeedNtHash: k should be equal to seed string lengths");
    get_blocks(seeds, blocks, monomers);
  }

  /**
   * Construct an ntHash object for spaced seeds.
   * @param seq String of the sequence to be hashed
   * @param seeds Vector of spaced seed patterns as strings (1s as cares, 0s
   * as don't cares, must be of size \p k)
   * @param num_hashes_per_seed Number of hashes to generate per seed
   * @param k K-mer size
   * @param pos Position in seq to start hashing from
   */
  SeedNtHash(const std::string& seq,
             const std::vector<std::string>& seeds,
             hashing_internals::NUM_HASHES_TYPE num_hashes_per_seed,
             hashing_internals::K_TYPE k,
             size_t pos = 0)
    : SeedNtHash(seq.data(), seq.size(), seeds, num_hashes_per_seed, k, pos)
  {
  }

  /**
   * Construct an ntHash object for spaced seeds.
   * @param seq C-string of the sequence to be hashed
   * @param seq_len Length of the sequence
   * @param seeds Vector of parsed spaced seed patterns (vectors of don't care
   * positions)
   * @param num_hashes_per_seed Number of hashes to generate per seed
   * @param k K-mer size
   * @param pos Position in seq to start hashing from
   */
  SeedNtHash(const char* seq,
             size_t seq_len,
             const std::vector<std::vector<unsigned>>& seeds,
             hashing_internals::NUM_HASHES_TYPE num_hashes_per_seed,
             hashing_internals::K_TYPE k,
             size_t pos = 0)
    : seq(seq)
    , seq_len(seq_len)
    , num_hashes_per_seed(num_hashes_per_seed)
    , k(k)
    , pos(pos)
    , initialized(false)
    , fwd_hash_nomonos(new uint64_t[seeds.size()])
    , rev_hash_nomonos(new uint64_t[seeds.size()])
    , fwd_hash(new uint64_t[seeds.size()])
    , rev_hash(new uint64_t[seeds.size()])
    , hash_arr(new uint64_t[num_hashes_per_seed * seeds.size()])
  {
    parsed_seeds_to_blocks(seeds, k, blocks, monomers);
  }

  /**
   * Construct an ntHash object for spaced seeds.
   * @param seq String of the sequence to be hashed
   * @param seeds Vector of parsed spaced seed patterns (vectors of don't care
   * positions)
   * @param num_hashes_per_seed Number of hashes to generate per seed
   * @param k K-mer size
   * @param pos Position in seq to start hashing from
   */
  SeedNtHash(const std::string& seq,
             const std::vector<std::vector<unsigned>>& seeds,
             hashing_internals::NUM_HASHES_TYPE num_hashes_per_seed,
             hashing_internals::K_TYPE k,
             size_t pos = 0)
    : SeedNtHash(seq.data(), seq.size(), seeds, num_hashes_per_seed, k, pos)
  {
  }

  SeedNtHash(const SeedNtHash& obj)
    : seq(obj.seq)
    , seq_len(obj.seq_len)
    , num_hashes_per_seed(obj.num_hashes_per_seed)
    , k(obj.k)
    , pos(obj.pos)
    , initialized(obj.initialized)
    , blocks(obj.blocks)
    , monomers(obj.monomers)
    , fwd_hash_nomonos(new uint64_t[obj.blocks.size()])
    , rev_hash_nomonos(new uint64_t[obj.blocks.size()])
    , fwd_hash(new uint64_t[obj.blocks.size()])
    , rev_hash(new uint64_t[obj.blocks.size()])
    , hash_arr(new uint64_t[obj.num_hashes_per_seed * obj.blocks.size()])
  {
    std::memcpy(fwd_hash_nomonos.get(),
                obj.fwd_hash_nomonos.get(),
                obj.blocks.size() * sizeof(uint64_t));
    std::memcpy(rev_hash_nomonos.get(),
                obj.rev_hash_nomonos.get(),
                obj.blocks.size() * sizeof(uint64_t));
    std::memcpy(
      fwd_hash.get(), obj.fwd_hash.get(), obj.blocks.size() * sizeof(uint64_t));
    std::memcpy(
      rev_hash.get(), obj.rev_hash.get(), obj.blocks.size() * sizeof(uint64_t));
    std::memcpy(hash_arr.get(),
                obj.hash_arr.get(),
                obj.num_hashes_per_seed * obj.blocks.size() * sizeof(uint64_t));
  }

  SeedNtHash(SeedNtHash&&) = default;

  /**
   * Reset iterator on a new sequence. Useful for re-using NtHash objects.
   * @param seq New sequence for hashing
   */
  void set_seq(const std::string& seq, size_t pos = 0)
  {
    this->seq = seq.data();
    this->seq_len = seq.size();
    this->pos = pos;
    this->initialized = false;
    this->fwd_hash.reset(new uint64_t[blocks.size()]);
    this->rev_hash.reset(new uint64_t[blocks.size()]);
    this->hash_arr.reset(new uint64_t[blocks.size() * num_hashes_per_seed]);
    check_error(this->seq_len < k,
                "SeedNtHash: sequence length (" +
                  std::to_string(this->seq_len) + ") is smaller than k (" +
                  std::to_string(k) + ")");
    check_error(pos > this->seq_len - k,
                "SeedNtHash: passed position (" + std::to_string(pos) +
                  ") is larger than sequence length (" +
                  std::to_string(this->seq_len) + ")");
  }

  /**
   * Calculate the next hash value. Refer to \ref NtHash::roll() for more
   * information.
   * @return \p true on success and \p false otherwise.
   */
  bool roll()
  {
    if (!initialized) {
      return init();
    }
    if (pos >= seq_len - k) {
      return false;
    }
    if (SEED_TAB[(unsigned char)seq[pos + k]] == SEED_N) {
      pos += k;
      return init();
    }
    ntmsm64(seq + pos,
            blocks,
            monomers,
            k,
            blocks.size(),
            num_hashes_per_seed,
            fwd_hash_nomonos.get(),
            rev_hash_nomonos.get(),
            fwd_hash.get(),
            rev_hash.get(),
            hash_arr.get());
    ++pos;
    return true;
  }

  /**
   * Like the roll() function, but advance backwards.
   * @return \p true on success and \p false otherwise.
   */
  bool roll_back()
  {
    if (!initialized) {
      return init();
    }
    if (pos == 0) {
      return false;
    }
    if (SEED_TAB[(unsigned char)seq[pos - 1]] == SEED_N && pos >= k) {
      pos -= k;
      return init();
    }
    if (SEED_TAB[(unsigned char)seq[pos - 1]] == SEED_N) {
      return false;
    }
    ntmsm64l(seq + pos - 1,
             blocks,
             monomers,
             k,
             blocks.size(),
             num_hashes_per_seed,
             fwd_hash_nomonos.get(),
             rev_hash_nomonos.get(),
             fwd_hash.get(),
             rev_hash.get(),
             hash_arr.get());
    --pos;
    return true;
  }

  /**
   * Peeks the hash values as if roll() was called. Refer to
   * \ref NtHash::peek() for more information.
   * @return \p true on success and \p false otherwise.
   */
  bool peek()
  {
    if (pos >= seq_len - k) {
      return false;
    }
    return peek(seq[pos + k]);
  }

  /**
   * Like peek(), but as if roll_back() was called.
   * @return \p true on success and \p false otherwise.
   */
  bool peek_back()
  {
    if (pos == 0) {
      return false;
    }
    return peek_back(seq[pos - 1]);
  }

  /**
   * Like peek(), but as if roll(char char_in) was called.
   * @return \p true on success and \p false otherwise.
   */
  bool peek(char char_in)
  {
    if (!initialized) {
      return init();
    }
    const std::unique_ptr<uint64_t[]> fwd_hash_nomonos_cpy(
      new uint64_t[blocks.size()]);
    const std::unique_ptr<uint64_t[]> rev_hash_nomonos_cpy(
      new uint64_t[blocks.size()]);
    const std::unique_ptr<uint64_t[]> fwd_hash_cpy(new uint64_t[blocks.size()]);
    const std::unique_ptr<uint64_t[]> rev_hash_cpy(new uint64_t[blocks.size()]);
    std::memcpy(fwd_hash_nomonos_cpy.get(),
                fwd_hash_nomonos.get(),
                blocks.size() * sizeof(uint64_t));
    std::memcpy(rev_hash_nomonos_cpy.get(),
                rev_hash_nomonos.get(),
                blocks.size() * sizeof(uint64_t));
    std::memcpy(
      fwd_hash_cpy.get(), fwd_hash.get(), blocks.size() * sizeof(uint64_t));
    std::memcpy(
      rev_hash_cpy.get(), rev_hash.get(), blocks.size() * sizeof(uint64_t));
    ntmsm64(seq + pos,
            char_in,
            blocks,
            monomers,
            k,
            blocks.size(),
            num_hashes_per_seed,
            fwd_hash_nomonos_cpy.get(),
            rev_hash_nomonos_cpy.get(),
            fwd_hash_cpy.get(),
            rev_hash_cpy.get(),
            hash_arr.get());
    return true;
  }

  /**
   * Like peek(), but as if roll_back(char char_in) was called.
   * @return \p true on success and \p false otherwise.
   */
  bool peek_back(char char_in)
  {
    if (!initialized) {
      return init();
    }
    const std::unique_ptr<uint64_t[]> fwd_hash_nomonos_cpy(
      new uint64_t[blocks.size()]);
    const std::unique_ptr<uint64_t[]> rev_hash_nomonos_cpy(
      new uint64_t[blocks.size()]);
    const std::unique_ptr<uint64_t[]> fwd_hash_cpy(new uint64_t[blocks.size()]);
    const std::unique_ptr<uint64_t[]> rev_hash_cpy(new uint64_t[blocks.size()]);
    std::memcpy(fwd_hash_nomonos_cpy.get(),
                fwd_hash_nomonos.get(),
                blocks.size() * sizeof(uint64_t));
    std::memcpy(rev_hash_nomonos_cpy.get(),
                rev_hash_nomonos.get(),
                blocks.size() * sizeof(uint64_t));
    std::memcpy(
      fwd_hash_cpy.get(), fwd_hash.get(), blocks.size() * sizeof(uint64_t));
    std::memcpy(
      rev_hash_cpy.get(), rev_hash.get(), blocks.size() * sizeof(uint64_t));
    ntmsm64l(seq + pos - 1,
             char_in,
             blocks,
             monomers,
             k,
             blocks.size(),
             num_hashes_per_seed,
             fwd_hash_nomonos_cpy.get(),
             rev_hash_nomonos_cpy.get(),
             fwd_hash_cpy.get(),
             rev_hash_cpy.get(),
             hash_arr.get());
    return true;
  }

  /**
   * Get the array of current hash values (length = \p get_hash_num())
   * @return Pointer to the hash array
   */
  const uint64_t* hashes() const { return hash_arr.get(); }

  /**
   * Get the position of last hashed k-mer or the k-mer to be hashed if roll()
   * has never been called on this NtHash object.
   * @return Position of the most recently hashed k-mer's first base-pair
   */
  size_t get_pos() const { return pos; }

  /**
   * Get the length of the hash array.
   * @return Number of seeds times \p get_hash_num_per_seed()
   */
  unsigned get_hash_num() const { return num_hashes_per_seed * blocks.size(); }

  /**
   * Get the number of hashes generated per seed.
   * @return Number of hashes per seed
   */
  hashing_internals::NUM_HASHES_TYPE get_hash_num_per_seed() const
  {
    return num_hashes_per_seed;
  }

  /**
   * Get the length of the k-mers.
   * @return \p k
   */
  hashing_internals::K_TYPE get_k() const { return k; }

  /**
   * Get the hash values of the forward strand.
   * @return Array of forward hash value arrays for each seed
   */
  uint64_t* get_forward_hash() const { return fwd_hash.get(); }

  /**
   * Get the hash values of the reverse strand.
   * @return Array of reverse-complement hash value arrays for each seed
   */
  uint64_t* get_reverse_hash() const { return rev_hash.get(); }

private:
  const char* seq;
  size_t seq_len;
  hashing_internals::NUM_HASHES_TYPE num_hashes_per_seed;
  hashing_internals::K_TYPE k;
  size_t pos;
  bool initialized;
  std::vector<hashing_internals::SpacedSeedBlocks> blocks;
  std::vector<hashing_internals::SpacedSeedMonomers> monomers;
  std::unique_ptr<uint64_t[]> fwd_hash_nomonos;
  std::unique_ptr<uint64_t[]> rev_hash_nomonos;
  std::unique_ptr<uint64_t[]> fwd_hash;
  std::unique_ptr<uint64_t[]> rev_hash;
  std::unique_ptr<uint64_t[]> hash_arr;

  /**
   * Initialize the internal state of the iterator
   * @return \p true if successful, \p false otherwise
   */
  bool init()
  {
    unsigned pos_n = 0;
    while (pos < seq_len - k + 1 && !ntmsm64(seq + pos,
                                             blocks,
                                             monomers,
                                             k,
                                             blocks.size(),
                                             num_hashes_per_seed,
                                             fwd_hash_nomonos.get(),
                                             rev_hash_nomonos.get(),
                                             fwd_hash.get(),
                                             rev_hash.get(),
                                             pos_n,
                                             hash_arr.get())) {
      pos += pos_n + 1;
    }
    if (pos > seq_len - k) {
      return false;
    }
    initialized = true;
    return true;
  }
};

/**
 * Similar to the SeedNtHash class, but instead of rolling on a predefined
 * sequence, BlindSeedNtHash needs to be fed the new character on each roll.
 */
class BlindSeedNtHash
{

public:
  /**
   * Construct an ntHash object for hashing spaced seeds on-the-fly.
   * @param seq C-string of the data. Only the first \p k characters will be
   * used, starting from \p pos.
   * @param seeds Vector of parsed spaced seed patterns (vectors of don't care
   * positions)
   * @param num_hashes_per_seed Number of hashes to generate per seed
   * @param k K-mer size
   * @param pos Position in seq to start hashing from
   */
  BlindSeedNtHash(const char* seq,
                  const std::vector<std::string>& seeds,
                  hashing_internals::NUM_HASHES_TYPE num_hashes_per_seed,
                  hashing_internals::K_TYPE k,
                  long pos = 0)
    : seq(seq + pos, seq + pos + k)
    , num_hashes_per_seed(num_hashes_per_seed)
    , k(k)
    , pos(pos)
    , fwd_hash_nomonos(new uint64_t[seeds.size()])
    , rev_hash_nomonos(new uint64_t[seeds.size()])
    , fwd_hash(new uint64_t[seeds.size()])
    , rev_hash(new uint64_t[seeds.size()])
    , hash_arr(new uint64_t[num_hashes_per_seed * seeds.size()])
  {
    check_seeds(seeds, k);
    get_blocks(seeds, blocks, monomers);
    unsigned pos_n = 0;
    ntmsm64(seq + pos,
            blocks,
            monomers,
            k,
            blocks.size(),
            num_hashes_per_seed,
            fwd_hash_nomonos.get(),
            rev_hash_nomonos.get(),
            fwd_hash.get(),
            rev_hash.get(),
            pos_n,
            hash_arr.get());
  }

  BlindSeedNtHash(const BlindSeedNtHash& seed_nthash)
    : seq(seed_nthash.seq)
    , num_hashes_per_seed(seed_nthash.num_hashes_per_seed)
    , k(seed_nthash.k)
    , pos(seed_nthash.pos)
    , blocks(seed_nthash.blocks)
    , monomers(seed_nthash.monomers)
    , fwd_hash_nomonos(new uint64_t[seed_nthash.blocks.size()])
    , rev_hash_nomonos(new uint64_t[seed_nthash.blocks.size()])
    , fwd_hash(new uint64_t[seed_nthash.blocks.size()])
    , rev_hash(new uint64_t[seed_nthash.blocks.size()])
    , hash_arr(new uint64_t[num_hashes_per_seed * seed_nthash.blocks.size()])
  {
    std::memcpy(fwd_hash_nomonos.get(),
                seed_nthash.fwd_hash_nomonos.get(),
                seed_nthash.blocks.size() * sizeof(uint64_t));
    std::memcpy(rev_hash_nomonos.get(),
                seed_nthash.rev_hash_nomonos.get(),
                seed_nthash.blocks.size() * sizeof(uint64_t));
    std::memcpy(fwd_hash.get(),
                seed_nthash.fwd_hash.get(),
                seed_nthash.blocks.size() * sizeof(uint64_t));
    std::memcpy(rev_hash.get(),
                seed_nthash.rev_hash.get(),
                seed_nthash.blocks.size() * sizeof(uint64_t));
    std::memcpy(hash_arr.get(),
                seed_nthash.hash_arr.get(),
                num_hashes_per_seed * seed_nthash.blocks.size() *
                  sizeof(uint64_t));
  }

  BlindSeedNtHash(BlindSeedNtHash&&) = default;

  /**
   * Like the NtHash::roll() function, but instead of advancing in the
   * sequence BlindSeedNtHash object was constructed on, the provided character
   * \p char_in is used as the next base.
   */
  void roll(char char_in)
  {
    seq.push_back(char_in);
    ntmsm64(seq,
            blocks,
            monomers,
            k,
            blocks.size(),
            num_hashes_per_seed,
            fwd_hash_nomonos.get(),
            rev_hash_nomonos.get(),
            fwd_hash.get(),
            rev_hash.get(),
            hash_arr.get());
    seq.pop_front();
    ++pos;
  }

  /**
   * Like the roll(char char_in) function, but advance backwards.
   */
  void roll_back(char char_in)
  {
    seq.push_front(char_in);
    ntmsm64l(seq,
             blocks,
             monomers,
             k,
             blocks.size(),
             num_hashes_per_seed,
             fwd_hash_nomonos.get(),
             rev_hash_nomonos.get(),
             fwd_hash.get(),
             rev_hash.get(),
             hash_arr.get());
    seq.pop_back();
    --pos;
  }

  /**
   * Get the array of current hash values (length = \p get_hash_num())
   * @return Pointer to the hash array
   */
  const uint64_t* hashes() const { return hash_arr.get(); }

  /**
   * Get the position of last hashed k-mer or the k-mer to be hashed if roll()
   * has never been called on this NtHash object.
   * @return Position of the most recently hashed k-mer's first base-pair
   */
  long get_pos() const { return pos; }

  /**
   * Get the length of the hash array.
   * @return Number of seeds times \p get_hash_num_per_seed()
   */
  unsigned get_hash_num() const { return num_hashes_per_seed * blocks.size(); }

  /**
   * Get the number of hashes generated per seed.
   * @return Number of hashes per seed
   */
  hashing_internals::NUM_HASHES_TYPE get_hash_num_per_seed() const
  {
    return num_hashes_per_seed;
  }

  /**
   * Get the length of the k-mers.
   * @return \p k
   */
  hashing_internals::K_TYPE get_k() const { return k; }

  /**
   * Get the hash values of the forward strand.
   * @return Array of forward hash value arrays for each seed
   */
  uint64_t* get_forward_hash() const { return fwd_hash.get(); }

  /**
   * Get the hash values of the reverse strand.
   * @return Array of reverse-complement hash value arrays for each seed
   */
  uint64_t* get_reverse_hash() const { return rev_hash.get(); }

private:
  std::deque<char> seq;
  hashing_internals::NUM_HASHES_TYPE num_hashes_per_seed;
  hashing_internals::K_TYPE k;
  long pos;
  std::vector<hashing_internals::SpacedSeedBlocks> blocks;
  std::vector<hashing_internals::SpacedSeedMonomers> monomers;
  std::unique_ptr<uint64_t[]> fwd_hash_nomonos;
  std::unique_ptr<uint64_t[]> rev_hash_nomonos;
  std::unique_ptr<uint64_t[]> fwd_hash;
  std::unique_ptr<uint64_t[]> rev_hash;
  std::unique_ptr<uint64_t[]> hash_arr;
};

} // namespace btllib
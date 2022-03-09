/*
 * nthash_lowlevel.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */

#ifndef BTLLIB_NTHASH_LOWLEVEL_HPP
#define BTLLIB_NTHASH_LOWLEVEL_HPP

#include "nthash_consts.hpp"
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace btllib {

// define a data structure for spaced seeds
// TODO: Create a clearer structure for list-of-blocks
using SpacedSeed = std::vector<unsigned>;

inline uint64_t
srol(const uint64_t x)
{
  uint64_t m = ((x & 0x8000000000000000ULL) >> 30) | // NOLINT
               ((x & 0x100000000ULL) >> 32);         // NOLINT
  return ((x << 1) & 0xFFFFFFFDFFFFFFFFULL) | m;     // NOLINT
}

inline uint64_t
srol(const uint64_t x, const unsigned d)
{
  uint64_t v = (x << d) | (x >> (64 - d));
  uint64_t y = (v ^ (v >> 33)) &                                   // NOLINT
               (std::numeric_limits<uint64_t>::max() >> (64 - d)); // NOLINT
  return v ^ (y | (y << 33));                                      // NOLINT
}

inline uint64_t
sror(const uint64_t x)
{
  uint64_t m = ((x & 0x200000000ULL) << 30) | ((x & 1ULL) << 32); // NOLINT
  return ((x >> 1) & 0xFFFFFFFEFFFFFFFFULL) | m;                  // NOLINT
}

// forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
inline uint64_t
ntf64(const char* kmer_seq, const unsigned k)
{
  uint64_t h_val = 0;
  for (unsigned i = 0; i < k / 4; i++) {
    h_val = srol(h_val, 4);
    uint8_t curr_offset = 4 * i;
    uint8_t tetramer_loc =
      64 * CONVERT_TAB[(unsigned char)kmer_seq[curr_offset]] +     // NOLINT
      16 * CONVERT_TAB[(unsigned char)kmer_seq[curr_offset + 1]] + // NOLINT
      4 * CONVERT_TAB[(unsigned char)kmer_seq[curr_offset + 2]] +
      CONVERT_TAB[(unsigned char)kmer_seq[curr_offset + 3]];
    h_val ^= TETRAMER_TAB[tetramer_loc];
  }
  unsigned remainder = k % 4;
  h_val = srol(h_val, remainder);
  if (remainder == 3) {
    uint8_t trimer_loc =
      16 * CONVERT_TAB[(unsigned char)kmer_seq[k - 3]] + // NOLINT
      4 * CONVERT_TAB[(unsigned char)kmer_seq[k - 2]] +
      CONVERT_TAB[(unsigned char)kmer_seq[k - 1]];
    // cppcheck-suppress arrayIndexOutOfBounds
    h_val ^= TRIMER_TAB[trimer_loc];
  } else if (remainder == 2) {
    uint8_t dimer_loc = 4 * CONVERT_TAB[(unsigned char)kmer_seq[k - 2]] +
                        CONVERT_TAB[(unsigned char)kmer_seq[k - 1]];
    // cppcheck-suppress arrayIndexOutOfBounds
    h_val ^= DIMER_TAB[dimer_loc];
  } else if (remainder == 1) {
    h_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1]];
  }
  return h_val;
}

// reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
inline uint64_t
ntr64(const char* kmer_seq, const unsigned k)
{
  uint64_t h_val = 0;
  unsigned remainder = k % 4;
  if (remainder == 3) {
    uint8_t trimer_loc =
      16 * RC_CONVERT_TAB[(unsigned char)kmer_seq[k - 1]] + // NOLINT
      4 * RC_CONVERT_TAB[(unsigned char)kmer_seq[k - 2]] +
      RC_CONVERT_TAB[(unsigned char)kmer_seq[k - 3]];
    // cppcheck-suppress arrayIndexOutOfBounds
    h_val ^= TRIMER_TAB[trimer_loc];
  } else if (remainder == 2) {
    uint8_t dimer_loc = 4 * RC_CONVERT_TAB[(unsigned char)kmer_seq[k - 1]] +
                        RC_CONVERT_TAB[(unsigned char)kmer_seq[k - 2]];
    // cppcheck-suppress arrayIndexOutOfBounds
    h_val ^= DIMER_TAB[dimer_loc];
  } else if (remainder == 1) {
    h_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1] & CP_OFF];
  }
  for (unsigned i = 0; i < k / 4; i++) {
    h_val = srol(h_val, 4);
    uint8_t curr_offset = 4 * (k / 4 - i) - 1;
    uint8_t tetramer_loc =
      64 * RC_CONVERT_TAB[(unsigned char)kmer_seq[curr_offset]] +     // NOLINT
      16 * RC_CONVERT_TAB[(unsigned char)kmer_seq[curr_offset - 1]] + // NOLINT
      4 * RC_CONVERT_TAB[(unsigned char)kmer_seq[curr_offset - 2]] +
      RC_CONVERT_TAB[(unsigned char)kmer_seq[curr_offset - 3]];
    h_val ^= TETRAMER_TAB[tetramer_loc];
  }
  return h_val;
}

// forward-strand ntHash for sliding k-mers
inline uint64_t
ntf64(const uint64_t fh_val,
      const unsigned k,
      const unsigned char char_out,
      const unsigned char char_in)
{
  uint64_t h_val = srol(fh_val);
  h_val ^= SEED_TAB[char_in];
  h_val ^= MS_TAB(char_out, k);
  return h_val;
}

// reverse-complement ntHash for sliding k-mers
inline uint64_t
ntr64(const uint64_t rh_val,
      const unsigned k,
      const unsigned char char_out,
      const unsigned char char_in)
{
  uint64_t h_val = rh_val ^ MS_TAB(char_in & CP_OFF, k);
  h_val ^= SEED_TAB[char_out & CP_OFF];
  h_val = sror(h_val);
  return h_val;
}

// canonical ntBase
inline uint64_t
ntc64(const char* kmer_seq, const unsigned k)
{
  uint64_t fh_val = 0, rh_val = 0;
  fh_val = ntf64(kmer_seq, k);
  rh_val = ntr64(kmer_seq, k);
  return (rh_val < fh_val) ? rh_val : fh_val;
}

// canonical ntHash
inline uint64_t
ntc64(const char* kmer_seq,
      const unsigned k,
      uint64_t& fh_val,
      uint64_t& rh_val)
{
  fh_val = ntf64(kmer_seq, k);
  rh_val = ntr64(kmer_seq, k);
  return (rh_val < fh_val) ? rh_val : fh_val;
}

// canonical ntHash for sliding k-mers
inline uint64_t
ntc64(const unsigned char char_out,
      const unsigned char char_in,
      const unsigned k,
      uint64_t& fh_val,
      uint64_t& rh_val)
{
  fh_val = ntf64(fh_val, k, char_out, char_in);
  rh_val = ntr64(rh_val, k, char_out, char_in);
  return (rh_val < fh_val) ? rh_val : fh_val;
}

// forward-strand ntHash for sliding k-mers to the left
inline uint64_t
ntf64l(const uint64_t rh_val,
       const unsigned k,
       const unsigned char char_out,
       const unsigned char char_in)
{
  uint64_t h_val = rh_val ^ MS_TAB(char_in, k);
  h_val ^= SEED_TAB[char_out];
  h_val = sror(h_val);
  return h_val;
}

// reverse-complement ntHash for sliding k-mers to the left
inline uint64_t
ntr64l(const uint64_t fh_val,
       const unsigned k,
       const unsigned char char_out,
       const unsigned char char_in)
{
  uint64_t h_val = srol(fh_val);
  h_val ^= SEED_TAB[char_in & CP_OFF];
  h_val ^= MS_TAB(char_out & CP_OFF, k);
  return h_val;
}

// canonical ntHash for sliding k-mers to the left
inline uint64_t
ntc64l(const unsigned char char_out,
       const unsigned char char_in,
       const unsigned k,
       uint64_t& fh_val,
       uint64_t& rh_val)
{
  fh_val = ntf64l(fh_val, k, char_out, char_in);
  rh_val = ntr64l(rh_val, k, char_out, char_in);
  return (rh_val < fh_val) ? rh_val : fh_val;
}

// one extra hash for given base hash
inline uint64_t
nte64(const uint64_t h_val, const unsigned k, const unsigned i)
{
  uint64_t t_val = h_val;
  t_val *= (i ^ k * MULTISEED);
  t_val ^= t_val >> MULTISHIFT;
  return t_val;
}

// canonical multihash ntBase
inline void
ntmc64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       uint64_t* h_val)
{
  uint64_t b_val = 0, t_val = 0;
  b_val = ntc64(kmer_seq, k);
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * MULTISEED);
    t_val ^= t_val >> MULTISHIFT;
    h_val[i] = t_val;
  }
}

// canonical multihash ntHash
inline void
ntmc64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val)
{
  uint64_t b_val = 0, t_val = 0;
  b_val = ntc64(kmer_seq, k, fh_val, rh_val);
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * MULTISEED);
    t_val ^= t_val >> MULTISHIFT;
    h_val[i] = t_val;
  }
}

// canonical multihash ntHash for sliding k-mers
inline void
ntmc64(const unsigned char char_out,
       const unsigned char char_in,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val)
{
  uint64_t b_val = 0, t_val = 0;
  b_val = ntc64(char_out, char_in, k, fh_val, rh_val);
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * MULTISEED);
    t_val ^= t_val >> MULTISHIFT;
    h_val[i] = t_val;
  }
}

// canonical multihash ntHash for sliding k-mers
inline void
ntmc64l(const unsigned char char_out,
        const unsigned char char_in,
        const unsigned k,
        const unsigned m,
        uint64_t& fh_val,
        uint64_t& rh_val,
        uint64_t* h_val)
{
  uint64_t b_val = 0, t_val = 0;
  b_val = ntc64l(char_out, char_in, k, fh_val, rh_val);
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * MULTISEED);
    t_val ^= t_val >> MULTISHIFT;
    h_val[i] = t_val;
  }
}

/*
 * ignoring k-mers containing nonACGT using ntHash function
 */

// canonical ntBase
inline bool
ntc64(const char* kmer_seq, const unsigned k, uint64_t& h_val, unsigned& loc_n)
{
  h_val = 0;
  loc_n = 0;
  uint64_t fh_val = 0, rh_val = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (SEED_TAB[(unsigned char)kmer_seq[i]] == SEED_N) {
      loc_n = i;
      return false;
    }
    fh_val = srol(fh_val);
    fh_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = srol(rh_val);
    rh_val ^= SEED_TAB[(unsigned char)kmer_seq[i] & CP_OFF];
  }
  h_val = (rh_val < fh_val) ? rh_val : fh_val;
  return true;
}

// canonical multihash ntBase
inline bool
ntmc64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       unsigned& loc_n,
       uint64_t* h_val)
{
  uint64_t b_val = 0, t_val = 0, fh_val = 0, rh_val = 0;
  loc_n = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (SEED_TAB[(unsigned char)kmer_seq[i]] == SEED_N) {
      loc_n = i;
      return false;
    }
    fh_val = srol(fh_val);
    fh_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = srol(rh_val);
    rh_val ^= SEED_TAB[(unsigned char)kmer_seq[i] & CP_OFF];
  }
  b_val = (rh_val < fh_val) ? rh_val : fh_val;
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * MULTISEED);
    t_val ^= t_val >> MULTISHIFT;
    h_val[i] = t_val;
  }
  return true;
}

// canonical ntHash
inline bool
ntc64(const char* kmer_seq,
      const unsigned k,
      uint64_t& fh_val,
      uint64_t& rh_val,
      uint64_t& h_val,
      unsigned& loc_n)
{
  h_val = fh_val = rh_val = 0;
  loc_n = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (SEED_TAB[(unsigned char)kmer_seq[i]] == SEED_N) {
      loc_n = i;
      return false;
    }
    fh_val = srol(fh_val);
    fh_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = srol(rh_val);
    rh_val ^= SEED_TAB[(unsigned char)kmer_seq[i] & CP_OFF];
  }
  h_val = (rh_val < fh_val) ? rh_val : fh_val;
  return true;
}

// canonical multihash ntHash
inline bool
ntmc64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       unsigned& loc_n,
       uint64_t* h_val)
{
  fh_val = rh_val = 0;
  uint64_t b_val = 0, t_val = 0;
  loc_n = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (SEED_TAB[(unsigned char)kmer_seq[i]] == SEED_N) {
      loc_n = i;
      return false;
    }
    fh_val = srol(fh_val);
    fh_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = srol(rh_val);
    rh_val ^= SEED_TAB[(unsigned char)kmer_seq[i] & CP_OFF];
  }
  b_val = (rh_val < fh_val) ? rh_val : fh_val;
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * MULTISEED);
    t_val ^= t_val >> MULTISHIFT;
    h_val[i] = t_val;
  }
  return true;
}

// strand-aware canonical multihash ntHash
inline bool
ntmc64(const char* kmer_seq,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       unsigned& loc_n,
       uint64_t* h_val,
       bool& h_stn)
{
  fh_val = rh_val = 0;
  uint64_t b_val = 0, t_val = 0;
  loc_n = 0;
  for (int i = int(k - 1); i >= 0; i--) {
    if (SEED_TAB[(unsigned char)kmer_seq[i]] == SEED_N) {
      loc_n = i;
      return false;
    }
    fh_val = srol(fh_val);
    fh_val ^= SEED_TAB[(unsigned char)kmer_seq[k - 1 - i]];

    rh_val = srol(rh_val);
    rh_val ^= SEED_TAB[(unsigned char)kmer_seq[i] & CP_OFF];
  }
  h_stn = rh_val < fh_val;
  b_val = h_stn ? rh_val : fh_val;
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * MULTISEED);
    t_val ^= t_val >> MULTISHIFT;
    h_val[i] = t_val;
  }
  return true;
}

// starnd-aware canonical multihash ntHash for sliding k-mers
inline void
ntmc64(const unsigned char char_out,
       const unsigned char char_in,
       const unsigned k,
       const unsigned m,
       uint64_t& fh_val,
       uint64_t& rh_val,
       uint64_t* h_val,
       bool& h_stn)
{
  uint64_t b_val = 0, t_val = 0;
  b_val = ntc64(char_out, char_in, k, fh_val, rh_val);
  h_stn = rh_val < fh_val;
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * MULTISEED);
    t_val ^= t_val >> MULTISHIFT;
    h_val[i] = t_val;
  }
}

// masking canonical ntHash using spaced seed pattern
inline uint64_t
mask_hash(uint64_t& fk_val,
          uint64_t& rk_val,
          const char* seed_seq,
          const char* kmer_seq,
          const unsigned k)
{
  uint64_t fs_val = fk_val, rs_val = rk_val;
  for (unsigned i = 0; i < k; i++) {
    if (seed_seq[i] != '1') {
      fs_val ^= MS_TAB((unsigned char)kmer_seq[i], k - 1 - i);
      rs_val ^= MS_TAB((unsigned char)kmer_seq[i] & CP_OFF, i);
    }
  }
  return (rs_val < fs_val) ? rs_val : fs_val;
}

// replacing canonical ntHash with a substitution
inline void
sub_hash(uint64_t fh_val,
         uint64_t rh_val,
         const char* kmer_seq,
         const std::vector<unsigned>& positions,
         const std::vector<unsigned char>& new_bases,
         const unsigned k,
         const unsigned m,
         uint64_t* h_val)
{
  uint64_t b_val = 0, t_val = 0;

  for (size_t i = 0; i < positions.size(); i++) {
    const auto pos = positions[i];
    const auto new_base = new_bases[i];

    fh_val ^= MS_TAB((unsigned char)kmer_seq[pos], k - 1 - pos);
    fh_val ^= MS_TAB(new_base, k - 1 - pos);

    rh_val ^= MS_TAB((unsigned char)kmer_seq[pos] & CP_OFF, pos);
    rh_val ^= MS_TAB(new_base & CP_OFF, pos);
  }

  b_val = rh_val < fh_val ? rh_val : fh_val;
  h_val[0] = b_val;
  for (unsigned i = 1; i < m; i++) {
    t_val = b_val * (i ^ k * MULTISEED);
    t_val ^= t_val >> MULTISHIFT;
    h_val[i] = t_val;
  }
}

// Multi spaced seed ntHash with multiple hashes per seed
inline bool
ntmsm64(const char* kmer_seq,
        const std::vector<SpacedSeed>& seed_seq,
        const std::vector<std::vector<unsigned>>& monomers,
        const unsigned k,
        const unsigned m,
        const unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        unsigned& loc_n,
        uint64_t* h_val)
{
  uint64_t fh_seed, rh_seed;
  unsigned i_base, block_start, block_end;
  const SpacedSeed* seed = nullptr;
  for (unsigned i_seed = 0; i_seed < m; i_seed++) {
    seed = &seed_seq[i_seed];
    fh_seed = 0;
    rh_seed = 0;
    for (int i_block = 0; i_block < (int)seed->size() - 1; i_block += 2) {
      // cppcheck-suppress arrayIndexOutOfBounds
      // cppcheck-suppress stlOutOfBounds
      block_start = seed->at(i_block);
      // cppcheck-suppress arrayIndexOutOfBounds
      // cppcheck-suppress stlOutOfBounds
      block_end = seed->at(i_block + 1);
      for (unsigned pos = block_start; pos < block_end; pos++) {
        if (kmer_seq[pos] == SEED_N) {
          loc_n = pos;
          return false;
        }
        fh_seed ^= MS_TAB((unsigned char)kmer_seq[pos], k - 1 - pos);
        rh_seed ^= MS_TAB((unsigned char)kmer_seq[pos] & CP_OFF, pos);
      }
    }
    fh_nomonos[i_seed] = fh_seed;
    rh_nomonos[i_seed] = rh_seed;
    for (unsigned pos : monomers[i_seed]) {
      fh_seed ^= MS_TAB((unsigned char)kmer_seq[pos], k - 1 - pos);
      rh_seed ^= MS_TAB((unsigned char)kmer_seq[pos] & CP_OFF, pos);
    }
    fh_val[i_seed] = fh_seed;
    rh_val[i_seed] = rh_seed;
    i_base = i_seed * m2;
    h_val[i_base] = fh_seed < rh_seed ? fh_seed : rh_seed;
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {
      h_val[i_base + i_hash] = h_val[i_base] * (i_hash ^ k * MULTISEED);
      h_val[i_base + i_hash] ^= h_val[i_base + i_hash] >> MULTISHIFT;
    }
  }
  return true;
}

#define NTMSM64(CHAR_IN_HANDLING)                                              \
  unsigned char char_out, char_in;                                             \
  uint64_t fh_seed, rh_seed;                                                   \
  unsigned i_out, i_in, i_base;                                                \
  const SpacedSeed* seed = nullptr;                                            \
  for (unsigned i_seed = 0; i_seed < m; i_seed++) {                            \
    seed = &seed_seq[i_seed];                                                  \
    fh_seed = srol(fh_nomonos[i_seed]);                                        \
    rh_seed = rh_nomonos[i_seed];                                              \
    for (int i_block = 0; i_block < (int)seed->size() - 1; i_block += 2) {     \
      /* cppcheck-suppress arrayIndexOutOfBounds */                            \
      /* cppcheck-suppress stlOutOfBounds */                                   \
      i_out = seed->at(i_block);                                               \
      /* cppcheck-suppress arrayIndexOutOfBounds */                            \
      /* cppcheck-suppress stlOutOfBounds */                                   \
      i_in = seed->at(i_block + 1);                                            \
      char_out = (unsigned char)kmer_seq[i_out];                               \
      CHAR_IN_HANDLING                                                         \
      fh_seed ^= MS_TAB(char_out, k - i_out);                                  \
      fh_seed ^= MS_TAB(char_in, k - i_in);                                    \
      rh_seed ^= MS_TAB(char_out & CP_OFF, i_out);                             \
      rh_seed ^= MS_TAB(char_in & CP_OFF, i_in);                               \
    }                                                                          \
    rh_seed = sror(rh_seed);                                                   \
    fh_nomonos[i_seed] = fh_seed;                                              \
    rh_nomonos[i_seed] = rh_seed;                                              \
    for (unsigned pos : monomers[i_seed]) {                                    \
      fh_seed ^= MS_TAB((unsigned char)kmer_seq[pos + 1], k - 1 - pos);        \
      rh_seed ^= MS_TAB((unsigned char)kmer_seq[pos + 1] & CP_OFF, pos);       \
    }                                                                          \
    fh_val[i_seed] = fh_seed;                                                  \
    rh_val[i_seed] = rh_seed;                                                  \
    i_base = i_seed * m2;                                                      \
    h_val[i_base] = fh_seed < rh_seed ? fh_seed : rh_seed;                     \
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {                         \
      h_val[i_base + i_hash] = h_val[i_base] * (i_hash ^ k * MULTISEED);       \
      h_val[i_base + i_hash] ^= h_val[i_base + i_hash] >> MULTISHIFT;          \
    }                                                                          \
  }

#define NTMSM64L(CHAR_IN_HANDLING)                                             \
  unsigned char char_in;                                                       \
  unsigned char char_out;                                                      \
  for (unsigned i_seed = 0; i_seed < m; i_seed++) {                            \
    const SpacedSeed& seed = seed_seq[i_seed];                                 \
    uint64_t rh_seed = srol(rh_val[i_seed]);                                   \
    uint64_t fh_seed = fh_val[i_seed];                                         \
    for (unsigned i_block = 0; i_block < seed.size() - 1; i_block += 2) {      \
      /* cppcheck-suppress arrayIndexOutOfBounds */                            \
      /* cppcheck-suppress stlOutOfBounds */                                   \
      const unsigned i_in = seed[i_block];                                     \
      /* cppcheck-suppress arrayIndexOutOfBounds */                            \
      /* cppcheck-suppress stlOutOfBounds */                                   \
      const unsigned i_out = seed[i_block + 1];                                \
      char_out = (unsigned char)kmer_seq[i_out];                               \
      CHAR_IN_HANDLING                                                         \
      fh_seed ^= MS_TAB(char_out, k - i_out);                                  \
      fh_seed ^= MS_TAB(char_in, k - i_in);                                    \
      rh_seed ^= MS_TAB(char_out & CP_OFF, i_out);                             \
      rh_seed ^= MS_TAB(char_in & CP_OFF, i_in);                               \
    }                                                                          \
    fh_seed = sror(fh_seed);                                                   \
    fh_val[i_seed] = fh_seed;                                                  \
    rh_val[i_seed] = rh_seed;                                                  \
    unsigned i_base = i_seed * m2;                                             \
    h_val[i_base] = fh_seed < rh_seed ? fh_seed : rh_seed;                     \
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {                         \
      h_val[i_base + i_hash] = h_val[i_base] * (i_hash ^ k * MULTISEED);       \
      h_val[i_base + i_hash] ^= h_val[i_base + i_hash] >> MULTISHIFT;          \
    }                                                                          \
  }

// Multi spaced seed ntHash for sliding k-mers with multiple hashes per seed
inline void
ntmsm64(const char* kmer_seq,
        const std::vector<SpacedSeed>& seed_seq,
        const std::vector<std::vector<unsigned>>& monomers,
        const unsigned k,
        const unsigned m,
        const unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        uint64_t* h_val)
{
  NTMSM64(char_in = (unsigned char)kmer_seq[i_in];)
}

// Multi spaced seed ntHash for sliding k-mers with multiple hashes per seed
inline void
ntmsm64l(const char* kmer_seq,
         const std::vector<SpacedSeed>& seed_seq,
         const unsigned k,
         const unsigned m,
         const unsigned m2,
         uint64_t* fh_val,
         uint64_t* rh_val,
         uint64_t* h_val)
{
  NTMSM64L(char_in = (unsigned char)kmer_seq[i_in];)
}

// Multi spaced seed ntHash for sliding k-mers with multiple hashes per seed
inline void
ntmsm64(const char* kmer_seq,
        const char in,
        const std::vector<SpacedSeed>& seed_seq,
        const std::vector<std::vector<unsigned>>& monomers,
        const unsigned k,
        const unsigned m,
        const unsigned m2,
        uint64_t* fh_nomonos,
        uint64_t* rh_nomonos,
        uint64_t* fh_val,
        uint64_t* rh_val,
        uint64_t* h_val)
{
  NTMSM64(
    if (i_in > k - 1) { char_in = in; } else {
      char_in = (unsigned char)kmer_seq[i_in];
    })
}

// Multi spaced seed ntHash for sliding k-mers with multiple hashes per seed
inline void
ntmsm64l(const char* kmer_seq,
         const char in,
         const std::vector<SpacedSeed>& seed_seq,
         const unsigned k,
         const unsigned m,
         const unsigned m2,
         uint64_t* fh_val,
         uint64_t* rh_val,
         uint64_t* h_val)
{
  NTMSM64L(
    if (i_in > k - 1) { char_in = in; } else {
      char_in = (unsigned char)kmer_seq[i_in];
    })
}

} // namespace btllib

#endif
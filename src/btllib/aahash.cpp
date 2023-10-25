#include <stdexcept>

#include "btllib/aahash.hpp"
#include <btllib/hashing_internals.hpp>

namespace {

using btllib::hashing_internals::LEVEL_X_AA_SEED_LEFT_31BITS_ROLL_TABLE;
using btllib::hashing_internals::LEVEL_X_AA_SEED_RIGHT_33BITS_ROLL_TABLE;
using btllib::hashing_internals::LEVEL_X_AA_SEED_TABLE;
using btllib::hashing_internals::srol;

inline uint64_t
base_hash(const char* kmer_seq, unsigned k, unsigned level = 1)
{
  uint64_t hash_value = 0;
  for (unsigned i = 0; i < k; i++) {
    hash_value = srol(hash_value);
    hash_value ^= LEVEL_X_AA_SEED_TABLE[level][(unsigned char)kmer_seq[i]];
  }
  return hash_value;
}

inline uint64_t
roll_hash(uint64_t hash_value,
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

inline uint64_t
modify_base_with_seed(uint64_t hash_value,
                      const btllib::SpacedSeed& seed,
                      const char* kmer_seq,
                      const unsigned k)
{
  for (unsigned i = 0; i < k; i++) {
    if (seed[i] != 1) {
      hash_value ^= AA_ROLL_TABLE((unsigned char)kmer_seq[i], 1, k - 1 - i);
      if (seed[i] != 0) {
        const unsigned level = seed[i];
        hash_value ^=
          AA_ROLL_TABLE((unsigned char)kmer_seq[i], level, k - 1 - i);
      }
    }
  }
  return hash_value;
}

} // namespace

namespace btllib {

using hashing_internals::AA_SEED__;
using hashing_internals::AA_SEED_TABLE;
using hashing_internals::extend_hashes;

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
  const uint64_t hash_value = base_hash(seq + pos, k, level);
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
    pos = std::numeric_limits<std::size_t>::max();
    return false;
  }
  if (AA_SEED_TABLE[(unsigned char)(seq[pos + k])] == AA_SEED__) {
    pos += k;
    return init();
  }
  const uint64_t hash_value =
    roll_hash(hashes_array[0], k, seq[pos], seq[pos + k], level);
  extend_hashes(hash_value, k, hash_num, hashes_array.get());
  ++pos;
  return true;
}

bool
SeedAAHash::roll()
{
  if (!aahash.roll()) {
    return false;
  }

  for (size_t i = 0; i < seeds.size(); ++i) {
    hashes_array[i * hash_num_per_seed] = modify_base_with_seed(
      aahash.hashes_array[0], seeds[i], aahash.seq + aahash.pos, aahash.k);
    extend_hashes(hashes_array[i * hash_num_per_seed],
                  aahash.k,
                  hash_num_per_seed,
                  hashes_array.get() + i * hash_num_per_seed);
  }

  return true;
}

bool
SeedAAHash::verify_seed()
{
  for (const auto& seed : seeds) {
    for (const auto& c : seed) {
      if (c != 0 && c != 1 && c != 2 && c != 3) {
        return false;
      }
    }
  }
  return true;
}

void
SeedAAHash::init()
{
  for (const auto& seed : seeds) {
    if (seed.size() != aahash.k) {
      throw std::runtime_error("Invalid seed. Seed length must be equal to k.");
    }
  }
  if (!verify_seed()) {
    throw std::runtime_error(
      "Invalid seed. Seed values must be 0, 1, 2, or 3.");
  }
}
} // namespace btllib

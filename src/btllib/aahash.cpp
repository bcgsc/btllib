#include <stdexcept>

#include "btllib/aahash.hpp"
#include "btllib/aahash_consts.hpp"
#include "btllib/nthash_lowlevel.hpp"

namespace btllib {
std::vector<SpacedSeed>
aa_parse_seeds(const std::vector<std::string>& seeds)
{
  std::vector<SpacedSeed> seed_vec;
  // convert vector of string to vector of vector of unsigned
  for (const auto& seed : seeds) {
    SpacedSeed seed_vec_tmp;
    for (const auto& c : seed) {
      // push back the unsigned value of the char by type conversion
      seed_vec_tmp.push_back((unsigned)(c - '0'));
    }
    seed_vec.push_back(seed_vec_tmp);
  }
  return seed_vec;
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
  const uint64_t hash_value = aahash_base(seq + pos, k, level);
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
    pos = std::numeric_limits<std::size_t>::max();
    return false;
  }
  if (AA_SEED_TABLE[(unsigned char)(seq[pos + k])] == AA_SEED__) {
    pos += k;
    return init();
  }
  const uint64_t hash_value =
    aahash_roll(hashes_array[0], k, seq[pos], seq[pos + k], level);
  nte64(hash_value, k, hash_num, hashes_array.get());
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
    hashes_array[i * hash_num_per_seed] = aa_modify_base_with_seed(
      aahash.hashes_array[0], seeds[i], aahash.seq + aahash.pos, aahash.k);
    nte64(hashes_array[i * hash_num_per_seed],
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

#pragma once

#include <array>
#include <cstdint>
#include <cstring>
#include <deque>
#include <memory>
#include <string>
#include <vector>

#include <btllib/hashing_internals.hpp>
#include <btllib/status.hpp>

namespace btllib::hashing_internals {

/**
 * Check the current k-mer for non ACGTU's
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return `true` if any of the first k characters is not an ACGTU, `false`
 * otherwise
 */
inline bool
is_invalid_kmer(const char* seq, unsigned k, size_t& pos_n)
{
  for (int i = (int)k - 1; i >= 0; i--) {
    if (SEED_TAB[(unsigned char)seq[i]] == SEED_N) {
      pos_n = i;
      return true;
    }
  }
  return false;
}

/**
 * Generate the forward-strand hash value of the first k-mer in the sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of k-mer_0
 */
inline uint64_t
base_forward_hash(const char* seq, unsigned k)
{
  uint64_t h_val = 0;
  for (unsigned i = 0; i < k - 3; i += 4) {
    h_val = srol(h_val, 4);
    uint8_t loc = 0;
    loc += 64 * CONVERT_TAB[(unsigned char)seq[i]];     // NOLINT
    loc += 16 * CONVERT_TAB[(unsigned char)seq[i + 1]]; // NOLINT
    loc += 4 * CONVERT_TAB[(unsigned char)seq[i + 2]];
    loc += CONVERT_TAB[(unsigned char)seq[i + 3]];
    h_val ^= TETRAMER_TAB[loc];
  }
  const unsigned remainder = k % 4;
  h_val = srol(h_val, remainder);
  if (remainder == 3) {
    uint8_t trimer_loc = 0;
    trimer_loc += 16 * CONVERT_TAB[(unsigned char)seq[k - 3]]; // NOLINT
    trimer_loc += 4 * CONVERT_TAB[(unsigned char)seq[k - 2]];
    trimer_loc += CONVERT_TAB[(unsigned char)seq[k - 1]];
    h_val ^= TRIMER_TAB[trimer_loc];
  } else if (remainder == 2) {
    uint8_t dimer_loc = 0;
    dimer_loc += 4 * CONVERT_TAB[(unsigned char)seq[k - 2]];
    dimer_loc += CONVERT_TAB[(unsigned char)seq[k - 1]];
    h_val ^= DIMER_TAB[dimer_loc];
  } else if (remainder == 1) {
    h_val ^= SEED_TAB[(unsigned char)seq[k - 1]];
  }
  return h_val;
}

/**
 * Perform a roll operation on the forward strand by removing char_out and
 * including char_in.
 * @param fh_val Previous hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Rolled forward hash value
 */
inline uint64_t
next_forward_hash(uint64_t fh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in)
{
  uint64_t h_val = srol(fh_val);
  h_val ^= SEED_TAB[char_in];
  h_val ^= srol_table(char_out, k);
  return h_val;
}

/**
 * Perform a roll back operation on the forward strand.
 * @param fh_val Previous hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Forward hash value rolled back
 */
inline uint64_t
prev_forward_hash(uint64_t fh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in)
{
  uint64_t h_val = fh_val ^ srol_table(char_in, k);
  h_val ^= SEED_TAB[char_out];
  h_val = sror(h_val);
  return h_val;
}

/**
 * Generate a hash value for the reverse-complement of the first k-mer in the
 * sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of the reverse-complement of k-mer_0
 */
inline uint64_t
base_reverse_hash(const char* seq, unsigned k)
{
  uint64_t h_val = 0;
  const unsigned remainder = k % 4;
  if (remainder == 3) {
    uint8_t trimer_loc = 0;
    trimer_loc += 16 * RC_CONVERT_TAB[(unsigned char)seq[k - 1]]; // NOLINT
    trimer_loc += 4 * RC_CONVERT_TAB[(unsigned char)seq[k - 2]];
    trimer_loc += RC_CONVERT_TAB[(unsigned char)seq[k - 3]];
    h_val ^= TRIMER_TAB[trimer_loc];
  } else if (remainder == 2) {
    uint8_t dimer_loc = 0;
    dimer_loc += 4 * RC_CONVERT_TAB[(unsigned char)seq[k - 1]];
    dimer_loc += RC_CONVERT_TAB[(unsigned char)seq[k - 2]];
    h_val ^= DIMER_TAB[dimer_loc];
  } else if (remainder == 1) {
    h_val ^= SEED_TAB[(unsigned char)seq[k - 1] & CP_OFF];
  }
  for (int i = (int)(k - remainder) - 1; i >= 3; i -= 4) {
    h_val = srol(h_val, 4);
    uint8_t loc = 0;
    loc += 64 * RC_CONVERT_TAB[(unsigned char)seq[i]];     // NOLINT
    loc += 16 * RC_CONVERT_TAB[(unsigned char)seq[i - 1]]; // NOLINT
    loc += 4 * RC_CONVERT_TAB[(unsigned char)seq[i - 2]];
    loc += RC_CONVERT_TAB[(unsigned char)seq[i - 3]];
    h_val ^= TETRAMER_TAB[loc];
  }
  return h_val;
}

/**
 * Perform a roll operation on the reverse-complement by removing char_out and
 * including char_in.
 * @param rh_val Previous reverse-complement hash value computed for the
 * sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Rolled hash value for the reverse-complement
 */
inline uint64_t
next_reverse_hash(uint64_t rh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in)
{
  uint64_t h_val = rh_val ^ srol_table(char_in & CP_OFF, k);
  h_val ^= SEED_TAB[char_out & CP_OFF];
  h_val = sror(h_val);
  return h_val;
}

/**
 * Perform a roll back operation on the reverse strand.
 * @param rh_val Previous hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Reverse hash value rolled back
 */
inline uint64_t
prev_reverse_hash(uint64_t rh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in)
{
  uint64_t h_val = srol(rh_val);
  h_val ^= SEED_TAB[char_in & CP_OFF];
  h_val ^= srol_table(char_out & CP_OFF, k);
  return h_val;
}

/**
 * Generate multiple new hash values for the input k-mer by substituting
 * multiple characters.
 *
 * @param fh_val Forward hash value of the k-mer.
 * @param rh_val Reverse hash value of the k-mer.
 * @param kmer_seq Array of characters representing the k-mer.
 * @param positions Indicies of the positions to be substituted.
 * @param new_bases Characters to be placed in the indicies indicated in
 * positions.
 * @param k k-mer size.
 * @param m Number of hashes per k-mer.
 * @param h_val Array of size m for storing the output hash values.
 */
inline void
sub_hash(uint64_t fh_val,
         uint64_t rh_val,
         const char* kmer_seq,
         const std::vector<unsigned>& positions,
         const std::vector<unsigned char>& new_bases,
         unsigned k,
         unsigned m,
         uint64_t* h_val)
{
  uint64_t b_val = 0;

  for (size_t i = 0; i < positions.size(); i++) {
    const auto pos = positions[i];
    const auto new_base = new_bases[i];

    fh_val ^= srol_table((unsigned char)kmer_seq[pos], k - 1 - pos);
    fh_val ^= srol_table(new_base, k - 1 - pos);

    rh_val ^= srol_table((unsigned char)kmer_seq[pos] & CP_OFF, pos);
    rh_val ^= srol_table(new_base & CP_OFF, pos);
  }

  b_val = canonical(fh_val, rh_val);
  extend_hashes(b_val, k, m, h_val);
}

} // namespace btllib::hashing_internals

namespace btllib {

using hashing_internals::base_forward_hash;
using hashing_internals::base_reverse_hash;
using hashing_internals::extend_hashes;
using hashing_internals::is_invalid_kmer;
using hashing_internals::next_forward_hash;
using hashing_internals::next_reverse_hash;
using hashing_internals::prev_forward_hash;
using hashing_internals::prev_reverse_hash;
using hashing_internals::SEED_N;
using hashing_internals::SEED_TAB;
using hashing_internals::sub_hash;

/**
 * Normal k-mer hashing.
 */
class NtHash
{

public:
  /**
   * Construct an ntHash object for k-mers.
   * @param seq C-string containing sequence data
   * @param seq_len Length of the sequence
   * @param num_hashes Number of hashes to generate per k-mer
   * @param k K-mer size
   * @param pos Position in the sequence to start hashing from
   */
  NtHash(const char* seq,
         size_t seq_len,
         hashing_internals::NUM_HASHES_TYPE num_hashes,
         hashing_internals::K_TYPE k,
         size_t pos = 0)
    : seq(seq)
    , seq_len(seq_len)
    , num_hashes(num_hashes)
    , k(k)
    , pos(pos)
    , initialized(false)
    , hash_arr(new uint64_t[num_hashes])
  {
    check_error(k == 0, "NtHash: k must be greater than 0");
    check_error(this->seq_len < k,
                "NtHash: sequence length (" + std::to_string(this->seq_len) +
                  ") is smaller than k (" + std::to_string(k) + ")");
    check_error(pos > this->seq_len - k,
                "NtHash: passed position (" + std::to_string(pos) +
                  ") is larger than sequence length (" +
                  std::to_string(this->seq_len) + ")");
  }

  /**
   * Construct an ntHash object for k-mers.
   * @param seq Sequence string
   * @param num_hashes Number of hashes to produce per k-mer
   * @param k K-mer size
   * @param pos Position in sequence to start hashing from
   */
  NtHash(const std::string& seq,
         hashing_internals::NUM_HASHES_TYPE num_hashes,
         hashing_internals::K_TYPE k,
         size_t pos = 0)
    : NtHash(seq.data(), seq.size(), num_hashes, k, pos)
  {
  }

  NtHash(const NtHash& obj)
    : seq(obj.seq)
    , seq_len(obj.seq_len)
    , num_hashes(obj.num_hashes)
    , k(obj.k)
    , pos(obj.pos)
    , initialized(obj.initialized)
    , fwd_hash(obj.fwd_hash)
    , rev_hash(obj.rev_hash)
    , hash_arr(new uint64_t[obj.num_hashes])
  {
    std::memcpy(
      hash_arr.get(), obj.hash_arr.get(), num_hashes * sizeof(uint64_t));
  }

  NtHash(NtHash&&) = default;

  /**
   * Calculate the hash values of current k-mer and advance to the next k-mer.
   * NtHash advances one nucleotide at a time until it finds a k-mer with valid
   * characters (ACGTU) and skips over those with invalid characters (non-ACGTU,
   * including N). This method must be called before hashes() is accessed, for
   * the first and every subsequent hashed kmer. get_pos() may be called at any
   * time to obtain the position of last hashed k-mer or the k-mer to be hashed
   * if roll() has never been called on this NtHash object. It is important to
   * note that the number of roll() calls is NOT necessarily equal to get_pos(),
   * if there are N's or invalid characters in the hashed sequence.
   * @return \p true on success and \p false otherwise
   */
  bool roll()
  {
    if (!initialized) {
      return init();
    }
    if (pos >= seq_len - k) {
      return false;
    }
    if (hashing_internals::SEED_TAB[(unsigned char)seq[pos + k]] ==
        hashing_internals::SEED_N) {
      pos += k;
      return init();
    }
    fwd_hash = next_forward_hash(fwd_hash, k, seq[pos], seq[pos + k]);
    rev_hash = next_reverse_hash(rev_hash, k, seq[pos], seq[pos + k]);
    extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
    ++pos;
    return true;
  }

  /**
   * Like the roll() function, but advance backwards.
   * @return \p true on success and \p false otherwise
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
    fwd_hash = prev_forward_hash(fwd_hash, k, seq[pos + k - 1], seq[pos - 1]);
    rev_hash = prev_reverse_hash(rev_hash, k, seq[pos + k - 1], seq[pos - 1]);
    extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
    --pos;
    return true;
  }

  /**
   * Peeks the hash values as if roll() was called (without advancing the
   * NtHash object. The peeked hash values can be obtained through the
   * hashes() method.
   * @return \p true on success and \p false otherwise
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
   * @return \p true on success and \p false otherwise
   */
  bool peek_back()
  {
    if (pos == 0) {
      return false;
    }
    return peek_back(seq[pos - 1]);
  }

  /**
   * Peeks the hash values as if roll() was called for char_in (without
   * advancing the NtHash object. The peeked hash values can be obtained through
   * the hashes() method.
   * @return \p true on success and \p false otherwise
   */
  bool peek(char char_in)
  {
    if (!initialized) {
      return init();
    }
    if (SEED_TAB[(unsigned char)char_in] == SEED_N) {
      return false;
    }
    const uint64_t fwd = next_forward_hash(fwd_hash, k, seq[pos], char_in);
    const uint64_t rev = next_reverse_hash(rev_hash, k, seq[pos], char_in);
    extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
    return true;
  }

  /**
   * Like peek(), but as if roll_back on char_in was called.
   * @return \p true on success and \p false otherwise
   */
  bool peek_back(char char_in)
  {
    if (!initialized) {
      return init();
    }
    if (SEED_TAB[(unsigned char)char_in] == SEED_N) {
      return false;
    }
    const unsigned char char_out = seq[pos + k - 1];
    const uint64_t fwd = prev_forward_hash(fwd_hash, k, char_out, char_in);
    const uint64_t rev = prev_reverse_hash(rev_hash, k, char_out, char_in);
    extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
    return true;
  }

  void sub(const std::vector<unsigned>& positions,
           const std::vector<unsigned char>& new_bases)
  {
    sub_hash(fwd_hash,
             rev_hash,
             seq + pos,
             positions,
             new_bases,
             get_k(),
             get_hash_num(),
             hash_arr.get());
  }

  /**
   * Get the array of current canonical hash values (length = \p get_hash_num())
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
   * Get the number of hashes generated per k-mer.
   * @return Number of hashes per k-mer
   */
  hashing_internals::NUM_HASHES_TYPE get_hash_num() const { return num_hashes; }

  /**
   * Get the length of the k-mers.
   * @return \p k
   */
  hashing_internals::K_TYPE get_k() const { return k; }

  /**
   * Get the hash value of the forward strand.
   * @return Forward hash value
   */
  uint64_t get_forward_hash() const { return fwd_hash; }

  /**
   * Get the hash value of the reverse strand.
   * @return Reverse-complement hash value
   */
  uint64_t get_reverse_hash() const { return rev_hash; }

private:
  const char* seq;
  const unsigned seq_len;
  hashing_internals::NUM_HASHES_TYPE num_hashes;
  hashing_internals::K_TYPE k;
  size_t pos;
  bool initialized;
  uint64_t fwd_hash = 0;
  uint64_t rev_hash = 0;
  std::unique_ptr<uint64_t[]> hash_arr;

  /**
   * Initialize the internal state of the iterator
   * @return \p true if successful, \p false otherwise
   */
  bool init()
  {
    size_t pos_n = 0;
    while (pos <= seq_len - k + 1 && is_invalid_kmer(seq + pos, k, pos_n)) {
      pos += pos_n + 1;
    }
    if (pos > seq_len - k) {
      return false;
    }
    fwd_hash = base_forward_hash(seq + pos, k);
    rev_hash = base_reverse_hash(seq + pos, k);
    extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
    initialized = true;
    return true;
  }
};

/**
 * Similar to the NtHash class, but instead of rolling on a predefined sequence,
 * BlindNtHash needs to be fed the new character on each roll. This is useful
 * when traversing an implicit de Bruijn Graph, as we need to query all bases
 * to know the possible extensions.
 */
class BlindNtHash
{

public:
  /**
   * Construct an ntHash object for hashing k-mers on-the-fly.
   * @param seq Sequence data. Only the first \p k characters will be
   * used, starting from \p pos.
   * @param hash_num Number of hashes to generate per k-mer
   * @param k K-mer size
   * @param pos Position in sequence to start hashing from
   */
  BlindNtHash(const std::string& seq,
              hashing_internals::NUM_HASHES_TYPE num_hashes,
              hashing_internals::K_TYPE k,
              long pos = 0)
    : seq(seq.data() + pos, seq.data() + pos + k)
    , num_hashes(num_hashes)
    , pos(pos)
    , fwd_hash(base_forward_hash(seq.data(), k))
    , rev_hash(base_reverse_hash(seq.data(), k))
    , hash_arr(new uint64_t[num_hashes])
  {
    check_error(k == 0, "BlindNtHash: k must be greater than 0");
    extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
  }

  BlindNtHash(const BlindNtHash& obj)
    : seq(obj.seq)
    , num_hashes(obj.num_hashes)
    , pos(obj.pos)
    , fwd_hash(obj.fwd_hash)
    , rev_hash(obj.rev_hash)
    , hash_arr(new uint64_t[obj.num_hashes])
  {
    std::memcpy(
      hash_arr.get(), obj.hash_arr.get(), num_hashes * sizeof(uint64_t));
  }

  BlindNtHash(BlindNtHash&&) = default;

  /**
   * Like the NtHash::roll() function, but instead of advancing in the
   * sequence BlindNtHash object was constructed on, the provided character
   * \p char_in is used as the next base. Useful if you want to query for
   * possible paths in an implicit de Bruijn graph graph.
   */
  void roll(char char_in)
  {
    fwd_hash = next_forward_hash(fwd_hash, seq.size(), seq.front(), char_in);
    rev_hash = next_reverse_hash(rev_hash, seq.size(), seq.front(), char_in);
    extend_hashes(fwd_hash, rev_hash, seq.size(), num_hashes, hash_arr.get());
    seq.pop_front();
    seq.push_back(char_in);
    ++pos;
  }

  /**
   * Like the roll(char char_in) function, but advance backwards.
   */
  void roll_back(char char_in)
  {
    fwd_hash = prev_forward_hash(fwd_hash, seq.size(), seq.back(), char_in);
    rev_hash = prev_reverse_hash(rev_hash, seq.size(), seq.back(), char_in);
    extend_hashes(fwd_hash, rev_hash, seq.size(), num_hashes, hash_arr.get());
    seq.pop_back();
    seq.push_front(char_in);
    --pos;
  }

  /**
   * Like NtHash::peek(), but as if roll(char char_in) was called.
   */
  void peek(char char_in)
  {
    const hashing_internals::K_TYPE k = seq.size();
    const uint64_t fwd = next_forward_hash(fwd_hash, k, seq.front(), char_in);
    const uint64_t rev = next_reverse_hash(rev_hash, k, seq.front(), char_in);
    extend_hashes(fwd, rev, seq.size(), num_hashes, hash_arr.get());
  }

  /**
   * Like peek(char char_in), but as if roll_back(char char_in) was called.
   */
  void peek_back(char char_in)
  {
    const hashing_internals::K_TYPE k = seq.size();
    const uint64_t fwd = prev_forward_hash(fwd_hash, k, seq.back(), char_in);
    const uint64_t rev = prev_reverse_hash(rev_hash, k, seq.back(), char_in);
    extend_hashes(fwd, rev, seq.size(), num_hashes, hash_arr.get());
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
   * Get the number of hashes generated per k-mer.
   * @return Number of hashes per k-mer
   */
  hashing_internals::NUM_HASHES_TYPE get_hash_num() const { return num_hashes; }

  /**
   * Get the length of the k-mers.
   * @return \p k
   */
  hashing_internals::K_TYPE get_k() const { return seq.size(); }

  /**
   * Get the hash value of the forward strand.
   * @return Forward hash value
   */
  uint64_t get_forward_hash() const { return fwd_hash; }

  /**
   * Get the hash value of the reverse strand.
   * @return Reverse-complement hash value
   */
  uint64_t get_reverse_hash() const { return rev_hash; }

private:
  std::deque<char> seq;
  hashing_internals::NUM_HASHES_TYPE num_hashes;
  long pos;
  uint64_t fwd_hash;
  uint64_t rev_hash;
  std::unique_ptr<uint64_t[]> hash_arr;
};

} // namespace btllib
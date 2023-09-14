#include "btllib/phred_nthash.hpp"

namespace btllib {
PhredNtHash::PhredNtHash(const char* seq,
                         size_t seq_len,
                         unsigned hash_num,
                         unsigned k,
                         size_t phred_min,
                         std::string_view quality_string,
                         size_t pos)
  : NtHash(seq, seq_len, hash_num, k, pos)
  , qual_seq(quality_string.data())
  , phred_min((unsigned char)(phred_min + (size_t)PHRED_OFFSET))
  , rmq(quality_string, seq_len)
{
}

PhredNtHash::PhredNtHash(const std::string& seq,
                         unsigned hash_num,
                         unsigned k,
                         size_t phred_min,
                         std::string_view quality_string,
                         size_t pos)
  : NtHash(seq, hash_num, k, pos)
  , qual_seq(quality_string.data())
  , phred_min((unsigned char)(phred_min + (size_t)PHRED_OFFSET))
  , rmq(quality_string, seq.length())
{
}

PhredNtHash::PhredNtHash(const char* seq,
                         size_t seq_len,
                         unsigned hash_num,
                         unsigned k,
                         size_t phred_min,
                         const char* quality_string,
                         size_t pos)
  : NtHash(seq, seq_len, hash_num, k, pos)
  , qual_seq(quality_string)
  , phred_min((unsigned char)(phred_min + (size_t)PHRED_OFFSET))
  , rmq(quality_string, seq_len)
{
}

PhredNtHash::PhredNtHash(const std::string& seq,
                         unsigned hash_num,
                         unsigned k,
                         size_t phred_min,
                         const char* quality_string,
                         size_t pos)
  : NtHash(seq, hash_num, k, pos)
  , qual_seq(quality_string)
  , phred_min((unsigned char)(phred_min + (size_t)PHRED_OFFSET))
  , rmq(quality_string, seq.length())
{
}

bool
PhredNtHash::roll()
{
  bool success = NtHash::roll();
  if (!success) {
    return false;
  }
  size_t curr_pos = NtHash::get_pos();
  size_t k = NtHash::get_k();
  size_t seq_len = NtHash::get_seq_len();
  size_t min_phred_idx = rmq.query(curr_pos, curr_pos + k - 1);
  auto min_phred = (unsigned char)qual_seq[min_phred_idx];
  while (min_phred < phred_min) {
    // check next kmer range does not exceed sequence length
    if (min_phred_idx + k >= seq_len) {
      return false;
    }
    curr_pos = min_phred_idx + 1;
    min_phred_idx = rmq.query(curr_pos, curr_pos + k - 1);
    min_phred = (size_t)qual_seq[min_phred_idx];
  }

  while (curr_pos != NtHash::get_pos()) {
    success = NtHash::roll();
  }

  return success;
}

bool
PhredNtHash::roll_back()
{
  bool success = NtHash::roll_back();
  if (!success) {
    return false;
  }
  size_t curr_pos = NtHash::get_pos();
  size_t k = NtHash::get_k();
  size_t min_phred_idx = rmq.query(curr_pos, curr_pos + k - 1);
  auto min_phred = (unsigned char)qual_seq[min_phred_idx];
  while (min_phred < phred_min) {
    // check next kmer range does not exceed sequence length
    if (min_phred_idx - k >= NtHash::get_seq_len()) {
      return false;
    }
    curr_pos = min_phred_idx - k;
    min_phred_idx = rmq.query(curr_pos, curr_pos + k - 1);
    min_phred = (size_t)qual_seq[min_phred_idx];
  }

  while (curr_pos != NtHash::get_pos()) {
    success = NtHash::roll_back();
  }

  return success;
}
} // namespace btllib
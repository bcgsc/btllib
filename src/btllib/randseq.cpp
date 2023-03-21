#include "btllib/randseq.hpp"
#include "btllib/status.hpp"

#include <algorithm>
#include <random>

namespace btllib {

RandSeq::RandSeq(SeqType seq_type, Masking masking)
{
  if (seq_type == SeqType::DNA) {
    chars = "ACGT";
  } else if (seq_type == SeqType::RNA) {
    chars = "ACGU";
  } else if (seq_type == SeqType::PROTEIN) {
    chars = "ACDEFGHIKLMNPQRSTVWY";
  }
  if (masking == Masking::SOFT) {
    std::string lowers = chars;
    std::transform(lowers.begin(),
                   lowers.end(),
                   lowers.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    chars += lowers;
  } else if (masking == Masking::HARD) {
    chars += seq_type == SeqType::PROTEIN ? 'X' : 'N';
  }
}

std::string
RandSeq::generate(size_t length)
{
  std::string seq;
  seq.reserve(length);
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<size_t> dist(0, chars.size() - 1);
  for (size_t i = 0; i < length; i++) {
    seq.append(std::string(1, chars[dist(rng)]));
  }
  return seq;
}

} // namespace btllib
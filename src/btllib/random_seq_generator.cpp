#include "btllib/random_seq_generator.hpp"
#include "btllib/status.hpp"

#include <algorithm>
#include <random>

namespace btllib {

RandomSequenceGenerator::RandomSequenceGenerator(SequenceType seq_type,
                                                 Masking masking)
{
  if (seq_type == DNA) {
    chars = "ACGT";
  } else if (seq_type == RNA) {
    chars = "ACGU";
  } else if (seq_type == PROTEIN) {
    chars = "ACDEFGHIKLMNPQRSTVWY";
  } else {
    check_error(true, CLASS_NAME + ": Invalid SequenceType");
  }
  if (masking == SOFT) {
    std::string lowers = chars;
    std::transform(lowers.begin(),
                   lowers.end(),
                   lowers.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    chars += lowers;
  } else if (masking == HARD) {
    chars += seq_type == PROTEIN ? 'X' : 'N';
  }
}

std::string
RandomSequenceGenerator::generate(size_t length)
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
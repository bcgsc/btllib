#ifndef BTLLIB_RANDSEQ_HPP
#define BTLLIB_RANDSEQ_HPP

#include <cstddef>
#include <string>

namespace btllib {

class RandSeq
{
public:
  enum class SeqType
  {
    DNA,
    RNA,
    PROTEIN
  };

  enum class Masking
  {
    NONE,
    SOFT,
    HARD
  };

  /**
   * Construct a random sequence generator object.
   *
   * @param type Sequence type (DNA, RNA, or protein)
   * @param masking If set to SOFT, lower-case values will also be generated. If
   * HARD, the sequences will include N/X positions.
   */
  RandSeq(SeqType type, Masking masking = Masking::NONE);

  /**
   * Generate a new random sequence.
   *
   * @param length Sequence length
   */
  std::string generate(size_t length);

private:
  std::string chars;
};

} // namespace btllib

#endif
#ifndef BTLLIB_RANDOM_SEQ_HPP
#define BTLLIB_RANDOM_SEQ_HPP

#include <cstddef>
#include <string>

namespace btllib {

class RandomSequenceGenerator
{
public:
  enum SequenceType
  {
    DNA,
    RNA,
    PROTEIN
  };

  enum Masking
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
  RandomSequenceGenerator(SequenceType type, Masking masking = NONE);

  /**
   * Generate a new random sequence.
   *
   * @param length Sequence length
   */
  std::string generate(size_t length);

private:
  const std::string CLASS_NAME = "RandomSequenceGenerator";
  std::string chars;
};

} // namespace btllib

#endif
#ifndef BTLLIB_RANDOM_SEQ_HPP
#define BTLLIB_RANDOM_SEQ_HPP

#include <stddef.h>
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

  RandomSequenceGenerator(SequenceType type, Masking masking = NONE);

  std::string generate(size_t length);

private:
  const std::string CLASS_NAME = "RandomSequenceGenerator";
  std::string chars;
};

} // namespace btllib

#endif
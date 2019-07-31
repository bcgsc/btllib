#ifndef BTL_SEQUENCE_READER_HPP
#define BTL_SEQUENCE_READER_HPP

#include "check.hpp"

#include <cstring>
#include <fstream>
#include <iostream>

namespace btl {

class SequenceReader
{
public:
  enum Flags
  {
    /** Fold lower-case characters to upper-case. */
    FOLD_CASE = 0,
    NO_FOLD_CASE = 1,
    /** Convert to standard quality. */
    NO_CONVERT_QUALITY = 0,
    CONVERT_QUALITY = 2,
  };

  SequenceReader(const char* filepath, int flags = 0);

private:
  std::ifstream ifs;
  std::istream& is;
  int flags;
};

inline SequenceReader::SequenceReader(const char* filepath, int flags)
  : ifs(filepath)
  , is(strcmp(filepath, "-") == 0 ? std::cin : ifs)
  , flags(flags)
{
  if (strcmp(filepath, "-") != 0) {
    check_stream(ifs, filepath);
  }
}

} // namespace btl

#endif
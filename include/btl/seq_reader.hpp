#ifndef BTL_SEQ_READER_HPP
#define BTL_SEQ_READER_HPP

#include "seq.hpp"
#include "status.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace btl {

/** Read a FASTA, FASTQ, SAM, or GFA2 file. */
class SeqReader
{
public:
  enum Flag
  {
    /** Fold lower-case characters to upper-case. */
    FOLD_CASE = 0,
    NO_FOLD_CASE = 1,
    /** Trim masked (lower case) characters from the ends of
     * sequences. */
    NO_TRIM_MASKED = 0,
    TRIM_MASKED = 2
  };

  SeqReader(const char* input_path, int flags = 0);

  bool flagFoldCase() const { return bool(~flags & NO_FOLD_CASE); }
  bool flagTrimMasked() const { return bool(flags & TRIM_MASKED); }

  enum Format
  {
    UNKNOWN,
    FASTA,
    FASTQ,
    SAM,
    GFA2,
    INVALID
  };

  Format get_format() const { return format; }

  /** Return the next character of this stream. */
  int peek() { return is.peek(); }

  /** Read operator. */
  bool read();

  /** Get the last read sequence. Cannot be called multiple times per read. */
  std::string seq();

  /** Quality of the last read sequence. Cannot be called multiple times per read. */
  std::string qual();

  /** Interface for manipulators. */
  void manip(std::istream& (*f)(std::istream&));

private:
  const char* input_path;
  std::ifstream ifs; // If reading from a file
  std::istream& is;  // Ref to the input stream (file or stdin)
  unsigned flags;
  Format format; // Format of the input file

  static const std::streamsize DETERMINE_FORMAT_MAX_READ_CHARS = 4096;

  static bool is_fasta(const char* input, size_t n);
  static bool is_fastq(const char* input, size_t n);
  static bool is_sam(const char* input, size_t n);
  static bool is_gfa2(const char* input, size_t n);

  void determine_format();

  void read_fasta();
  void read_fastq();
  void read_sam();
  void read_gfa2();

  void (SeqReader::*read_impl)();

  std::string s;
  std::string q;
  std::string tmp;
};

inline SeqReader::SeqReader(const char* input_path, int flags)
  : input_path(input_path)
  , ifs(input_path)
  , is(strcmp(input_path, "-") == 0 ? std::cin : ifs)
  , flags(flags)
  , format(Format::UNKNOWN)
  , read_impl(nullptr)
{
  if (strcmp(input_path, "-") != 0) {
    check_stream(ifs, input_path);
  }
  check_warning(is.eof(), std::string(input_path) + " is empty.");
  determine_format();
}

inline bool
SeqReader::is_fasta(const char* input, size_t n)
{
  size_t current = 0;
  unsigned char c;
  enum State
  {
    IN_HEADER_1,
    IN_HEADER_2,
    IN_SEQ
  };
  State state = IN_HEADER_1;
  while (current < n) {
    c = input[current];
    switch (state) {
      case IN_HEADER_1:
        if (c == '>') {
          state = IN_HEADER_2;
        } else {
          return false;
        }
        break;
      case IN_HEADER_2:
        if (c == '\n') {
          state = IN_SEQ;
        }
        break;
      case IN_SEQ:
        if (c == '\n') {
          state = IN_HEADER_1;
        } else if (!bool(COMPLEMENTS[c])) {
          return false;
        }
        break;
    }
    current++;
  }
  return true;
}

inline bool
SeqReader::is_fastq(const char* input, size_t n)
{
  size_t current = 0;
  unsigned char c;
  enum State
  {
    IN_HEADER_1,
    IN_HEADER_2,
    IN_SEQ,
    IN_PLUS_1,
    IN_PLUS_2,
    IN_QUAL
  };
  State state = IN_HEADER_1;
  while (current < n) {
    c = input[current];
    switch (state) {
      case IN_HEADER_1:
        if (c == '>') {
          state = IN_HEADER_2;
        } else {
          return false;
        }
        break;
      case IN_HEADER_2:
        if (c == '\n') {
          state = IN_SEQ;
        }
        break;
      case IN_SEQ:
        if (c == '\n') {
          state = IN_PLUS_1;
        } else if (!bool(COMPLEMENTS[c])) {
          return false;
        }
        break;
      case IN_PLUS_1:
        if (c == '+') {
          state = IN_PLUS_2;
        } else {
          return false;
        }
        break;
      case IN_PLUS_2:
        if (c == '\n') {
          state = IN_QUAL;
        }
        break;
      case IN_QUAL:
        if (c == '\n') {
          state = IN_HEADER_1;
        } else if (c < '!' || c > '~') {
          return false;
        }
        break;
    }
    current++;
  }
  return true;
}

inline bool
SeqReader::is_sam(const char* input, size_t n)
{
  enum Column
  {
    QNAME = 1,
    FLAG,
    RNAME,
    POS,
    MAPQ,
    CIGAR,
    RNEXT,
    PNEXT,
    TLEN,
    SEQ,
    QUAL
  };

  size_t current = 0;

  while (current < n && input[current] == '@') {
    while (current < n && input[current] != '\n') {
      current++;
    }
    current++;
  }

  int column = 1;
  unsigned char c;
  while (current < n) {
    c = input[current];
    if (c == '\n') {
      break;
    }
    if (c == '\t') {
      if (current > 0 && !bool(std::isspace(input[current - 1]))) {
        column++;
      } else {
        return false;
      }
    } else {
      switch (Column(column)) {
        case QNAME:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        case FLAG:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case RNAME:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        case POS:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case MAPQ:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case CIGAR:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        case RNEXT:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        case PNEXT:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case TLEN:
          if (!bool(std::isdigit(c))) {
            return false;
          }
          break;
        case SEQ:
          if (!bool(COMPLEMENTS[c])) {
            return false;
          }
          break;
        case QUAL:
          if (bool(std::isspace(c))) {
            return false;
          }
          break;
        default:
          break;
      }
    }
    current++;
  }

  return current >= n || column >= QUAL;
}

inline bool
SeqReader::is_gfa2(const char* input, size_t n) {
  const char specs[] = { 'H', 'S', 'F', 'E', 'G', 'O', 'U' };

  enum State
  {
    IN_ID,
    IN_ID_TAB,
    IN_REST,
    IN_IGNORED
  };

  auto is_a_spec = [&] (char c) {
    bool found = false;
    for (size_t i = 0; i < sizeof(specs) / sizeof(char); i++) {
      if (c == specs[i]) {
        found = true;
        break;
      }
    }
    return found;
  };

  State state = is_a_spec(input[0]) ? IN_ID : IN_IGNORED;
  bool has_id = false;
  size_t current = 0;
  unsigned char c;
  while (current < n) {
    c = input[current];
    switch (state) {
      case IN_ID:
        if (!is_a_spec(c)) { return false; }
        has_id = true;
        state = IN_ID_TAB;
        break;
      case IN_ID_TAB:
        if (c != '\t') { return false; }
        state = IN_REST;
        break;
      case IN_REST:
        if (c == '\n') {
          if (current + 1 < n) {
            state = is_a_spec(input[current + 1]) ? IN_ID : IN_IGNORED;
          }
        }
        break;
      case IN_IGNORED:
        if (c == '\n') {
          if (current + 1 < n) {
            state = is_a_spec(input[current + 1]) ? IN_ID : IN_IGNORED;
          }
        }
        break;
      default:
        break;
    }
    current++;
  }

  return has_id;
}

inline void
SeqReader::determine_format()
{
  std::ifstream ifs(input_path);
  std::istream& is(strcmp(input_path, "-") == 0 ? std::cin : ifs);

  char buffer[DETERMINE_FORMAT_MAX_READ_CHARS];
  is.read(buffer, DETERMINE_FORMAT_MAX_READ_CHARS);

  const auto n = is.gcount();
  if (is_fasta(buffer, n)) {
    format = Format::FASTA;
    read_impl = &SeqReader::read_fasta;
  } else if (is_fastq(buffer, n)) {
    format = Format::FASTQ;
    read_impl = &SeqReader::read_fastq;
  } else if (is_sam(buffer, n)) {
    format = Format::SAM;
    read_impl = &SeqReader::read_sam;
  } else if (is_gfa2(buffer, n)) {
    format = Format::GFA2;
    read_impl = &SeqReader::read_gfa2;
  } else {
    format = Format::INVALID;
    log_error(std::string(input_path) + " input file is in invalid format!");
    std::exit(EXIT_FAILURE);
  }
}

inline void
SeqReader::read_fasta()
{
  std::getline(is, tmp);
  std::getline(is, s);
}

inline void
SeqReader::read_fastq()
{
  std::getline(is, tmp);
  std::getline(is, s);
  std::getline(is, tmp);
  std::getline(is, q);
}

inline void
SeqReader::read_sam()
{
  enum Column
  {
    QNAME = 1,
    FLAG,
    RNAME,
    POS,
    MAPQ,
    CIGAR,
    RNEXT,
    PNEXT,
    TLEN,
    SEQ,
    QUAL
  };
  for (;;) {
    std::getline(is, tmp);
    if (tmp.length() > 0 && tmp[0] != '@') {
      size_t pos = 0, pos2 = 0, pos3 = 0;
      for (int i = 0; i < int(SEQ) - 1; i++) {
        pos = tmp.find('\t', pos + 1);
      }
      pos2 = tmp.find('\t', pos + 1);
      pos3 = tmp.find('\t', pos2 + 1);
      if (pos3 == std::string::npos) {
        pos3 = tmp.length();
      }

      s = tmp.substr(pos + 1, pos2 - pos - 1);
      q = tmp.substr(pos2 + 1, pos3 - pos2 - 1);

      break;
    }
    if (is.eof()) {
      break;
    }
  }
}

inline void
SeqReader::read_gfa2() {
  enum Column
  {
    S = 1, ID, LEN, SEQ
  };
  for (;;) {
    std::getline(is, tmp);
    if (tmp.length() > 0 && tmp[0] == 'S') {
      size_t pos = 0, pos2 = 0;
      for (int i = 0; i < int(SEQ) - 1; i++) {
        pos = tmp.find('\t', pos + 1);
      }
      pos2 = tmp.find('\t', pos + 1);
      if (pos2 == std::string::npos) {
        pos2 = tmp.length();
      }

      s = tmp.substr(pos + 1, pos2 - pos - 1);

      break;
    }
    if (is.eof()) {
      break;
    }
  }
}

inline bool
SeqReader::read()
{
  if (is.good() && is.peek() != std::istream::traits_type::eof()) {
    (this->*read_impl)();
    if (s.empty()) { return false; }
    if (flagTrimMasked()) {
      const auto len = s.length();
      size_t trim_start = 0, trim_end = s.length();
      while (trim_start <= len && bool(islower(s[trim_start]))) {
        trim_start++;
      }
      while (trim_end > 0 && bool(islower(s[trim_end - 1]))) {
        trim_end--;
      }
      s.erase(trim_end);
      s.erase(0, trim_start);
      if (!q.empty()) {
        q.erase(trim_end);
        q.erase(0, trim_start);
      }
    }
    if (flagFoldCase()) {
      std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    }
    return true;
  }
  return false;
}

inline std::string
SeqReader::seq()
{
  return std::move(s);
}

inline std::string
SeqReader::qual()
{
  return std::move(q);
}

inline void
SeqReader::manip(std::istream& (*f)(std::istream&))
{
  f(is);
}

} // namespace btl

#endif
#ifndef BTLLIB_SEQ_READER_HPP
#define BTLLIB_SEQ_READER_HPP

#include "data_saveload.hpp"
#include "seq.hpp"
#include "status.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <string>

namespace btllib {

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

  SeqReader(const std::string& source, int flags = 0);

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

  void close();

  Format get_format() const { return format; }

  /** Return the next character of this stream. */
  int peek();

  /** Read operator. */
  bool read();

  /** Last read sequence name. Cannot be called multiple times per read. */
  std::string name();

  /** Last read sequence comment. Cannot be called multiple times per read. */
  std::string comment();

  /** Last read sequence. Cannot be called multiple times per read. */
  std::string seq();

  /** Last read sequence quality. Cannot be called multiple times per
   * read. */
  std::string qual();

private:
  const std::string& source;
  DataSource input;
  bool closed = false;
  unsigned flags;
  Format format = UNKNOWN; // Format of the input file

  static const size_t RESERVE_SIZE_FOR_STRINGS = 1024;

  static const size_t BUFFER_SIZE = 2048;
  char buffer[BUFFER_SIZE];
  size_t buffer_start = 0;
  size_t buffer_end = 0;
  bool eof_newline_inserted = false;

  bool load_buffer();

  bool is_fasta();
  bool is_fastq();
  bool is_sam();
  bool is_gfa2();

  void determine_format();

  bool getline(std::string& line);
  char getsection(std::string& section);

  int read_stage = 0;
  bool read_fasta();
  bool read_fastq();
  bool read_sam();
  bool read_gfa2();

  bool (SeqReader::*read_impl)() = nullptr;

  std::string name_str;
  std::string comment_str;
  std::string seq_str;
  std::string qual_str;
  std::string tmp;
};

inline SeqReader::SeqReader(const std::string& source, int flags)
  : source(source)
  , input(data_load(source))
  , flags(flags)
{
  determine_format();
  name_str.reserve(RESERVE_SIZE_FOR_STRINGS);
  comment_str.reserve(RESERVE_SIZE_FOR_STRINGS);
  seq_str.reserve(RESERVE_SIZE_FOR_STRINGS);
  qual_str.reserve(RESERVE_SIZE_FOR_STRINGS);
  tmp.reserve(RESERVE_SIZE_FOR_STRINGS);
}

inline void
SeqReader::close()
{
  if (!closed) {
    input.close();
    closed = true;
  }
}

inline int
SeqReader::peek()
{
  int p = std::fgetc(input);
  std::ungetc(p, input);
  return p;
}

inline bool
SeqReader::load_buffer() {
  buffer_start = 0;
  while ((buffer_end = fread(buffer, 1, BUFFER_SIZE, input)) == 0 && !bool(std::feof(input))) {}
  if (bool(std::feof(input)) && !eof_newline_inserted) {
    if (buffer_end < BUFFER_SIZE) {
      if (buffer_end == 0 || (buffer_end > 0 && buffer[buffer_end - 1] != '\n')) {
        buffer[buffer_end++] = '\n';
      }
      eof_newline_inserted = true;
    } else if (buffer[BUFFER_SIZE - 1] == '\n') {
      eof_newline_inserted = true;
    }
    return true;
  }
  return bool(buffer_end);
}

inline bool
SeqReader::is_fasta()
{
  size_t current = buffer_start;
  unsigned char c;
  enum State
  {
    IN_HEADER_1,
    IN_HEADER_2,
    IN_SEQ
  };
  State state = IN_HEADER_1;
  while (current < buffer_end) {
    c = buffer[current];
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
SeqReader::is_fastq()
{
  size_t current = buffer_start;
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
  while (current < buffer_end) {
    c = buffer[current];
    switch (state) {
      case IN_HEADER_1:
        if (c == '@') {
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
SeqReader::is_sam()
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

  size_t current = buffer_start;

  while (current < buffer_end && buffer[current] == '@') {
    while (current < buffer_end && buffer[current] != '\n') {
      current++;
    }
    current++;
  }

  int column = 1;
  unsigned char c;
  while (current < buffer_end) {
    c = buffer[current];
    if (c == '\n') {
      break;
    }
    if (c == '\t') {
      if (current > 0 && !bool(std::isspace(buffer[current - 1]))) {
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

  return current >= buffer_end || column >= QUAL;
}

inline bool
SeqReader::is_gfa2()
{
  const unsigned char specs[] = { 'H', 'S', 'F', 'E', 'G', 'O', 'U' };

  enum State
  {
    IN_ID,
    IN_ID_TAB,
    IN_REST,
    IN_IGNORED
  };

  auto is_a_spec = [&](unsigned char c) {
    bool found = false;
    for (unsigned char spec : specs) {
      if (c == spec) {
        found = true;
        break;
      }
    }
    return found;
  };

  State state = is_a_spec(buffer[0]) ? IN_ID : IN_IGNORED;
  bool has_id = false;
  size_t current = buffer_start;
  unsigned char c;
  while (current < buffer_end) {
    c = buffer[current];
    switch (state) {
      case IN_ID:
        if (!is_a_spec(c)) {
          return false;
        }
        has_id = true;
        state = IN_ID_TAB;
        break;
      case IN_ID_TAB:
        if (c != '\t') {
          return false;
        }
        state = IN_REST;
        break;
      case IN_REST:
        if (c == '\n') {
          if (current + 1 < buffer_end) {
            state = is_a_spec(buffer[current + 1]) ? IN_ID : IN_IGNORED;
          }
        }
        break;
      case IN_IGNORED:
        if (c == '\n') {
          if (current + 1 < buffer_end) {
            state = is_a_spec(buffer[current + 1]) ? IN_ID : IN_IGNORED;
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
  load_buffer();
  check_warning(buffer_end - buffer_start == 1, std::string(source) + " is empty.");

  if (is_fasta()) {
    format = Format::FASTA;
    read_impl = &SeqReader::read_fasta;
  } else if (is_fastq()) {
    format = Format::FASTQ;
    read_impl = &SeqReader::read_fastq;
  } else if (is_sam()) {
    format = Format::SAM;
    read_impl = &SeqReader::read_sam;
  } else if (is_gfa2()) {
    format = Format::GFA2;
    read_impl = &SeqReader::read_gfa2;
  } else {
    format = Format::INVALID;
    log_error(std::string(source) + " input file is in invalid format!");
    std::exit(EXIT_FAILURE);
  }
}

inline bool
SeqReader::getline(std::string& line)
{
  char c;
  for (; buffer_start < buffer_end && (c = buffer[buffer_start]) != '\n'; ++buffer_start) {
    line += c;
  }
  if (c == '\n') {
    ++buffer_start;
    return true;
  }
  return false;
}

inline char
SeqReader::getsection(std::string& section)
{
  char c;
  for (; buffer_start < buffer_end && (c = buffer[buffer_start]) != '\n' && c != ' '; ++buffer_start) {
    section += c;
  }
  if (c == '\n' || c == ' ') {
    ++buffer_start;
    return c;
  }
  return char(0);
}

inline bool
SeqReader::read_fasta()
{
  switch (read_stage) {
    case 0: {
      if (!getline(tmp)) { return false; }
      auto pos = tmp.find(' ');
      name_str = tmp.substr(1, pos - 1);
      while (pos < tmp.size() && tmp[pos] == ' ') {
        pos++;
      }
      if (pos < tmp.size()) {
        comment_str = tmp.substr(pos);
      } else {
        comment_str = "";
      }
      ++read_stage;
      tmp.clear();
    }
    case 1: {
      if (!getline(tmp)) { return false; }
      seq_str = std::move(tmp);
      read_stage = 0;
      tmp.clear();
      return true;    
    }
  }
  return false;
}

inline bool
SeqReader::read_fastq()
{
  switch (read_stage) {
    case 0: {
      if (!getline(tmp)) { return false; }
      auto pos = tmp.find(' ');
      name_str = tmp.substr(1, pos - 1);
      while (pos < tmp.size() && tmp[pos] == ' ') {
        pos++;
      }
      if (pos < tmp.size()) {
        comment_str = tmp.substr(pos);
      } else {
        comment_str = "";
      }
      ++read_stage;
      tmp.clear();
    }
    case 1: {
      if (!getline(tmp)) { return false; }
      seq_str = std::move(tmp);
      ++read_stage;
      tmp.clear();
    }
    case 2: {
      if (!getline(tmp)) { return false; }
      ++read_stage;
      tmp.clear();
    }
    case 3: {
      if (!getline(tmp)) { return false; }
      qual_str = std::move(tmp);
      read_stage = 0;
      tmp.clear();
      return true;
    }
  }
  return false;
}

inline bool
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
    if (!getline(tmp)) { return false; }
    if (tmp.length() > 0 && tmp[0] != '@') {
      size_t pos = 0, pos2 = 0, pos3 = 0;
      pos2 = tmp.find('\t');
      name_str = tmp.substr(0, pos2);
      for (int i = 0; i < int(SEQ) - 1; i++) {
        pos = tmp.find('\t', pos + 1);
      }
      pos2 = tmp.find('\t', pos + 1);
      pos3 = tmp.find('\t', pos2 + 1);
      if (pos3 == std::string::npos) {
        pos3 = tmp.length();
      }

      seq_str = tmp.substr(pos + 1, pos2 - pos - 1);
      qual_str = tmp.substr(pos2 + 1, pos3 - pos2 - 1);

      tmp.clear();
      return true;
    }
    tmp.clear();
    if (buffer_start >= buffer_end) {
      return false;
    }
  }
}

inline bool
SeqReader::read_gfa2()
{
  enum Column
  {
    S = 1,
    ID,
    LEN,
    SEQ
  };
  for (;;) {
    if (!getline(tmp)) { return false; }
    if (tmp.length() > 0 && tmp[0] == 'S') {
      size_t pos = 0, pos2 = 0;
      pos2 = tmp.find('\t', 1);
      name_str = tmp.substr(1, pos2 - 1);
      for (int i = 0; i < int(SEQ) - 1; i++) {
        pos = tmp.find('\t', pos + 1);
      }
      pos2 = tmp.find('\t', pos + 1);
      if (pos2 == std::string::npos) {
        pos2 = tmp.length();
      }

      seq_str = tmp.substr(pos + 1, pos2 - pos - 1);

      tmp.clear();
      return true;
    }
    tmp.clear();
    if (buffer_start >= buffer_end) {
      return false;
    }
  }
}

inline bool
SeqReader::read()
{
  if (buffer_start < buffer_end) {
    while (!(this->*read_impl)()) {
      if (!load_buffer()) { return false; }
    }
    if (buffer_start >= buffer_end) {
      load_buffer();
    }

    if (seq_str.empty()) {
      return false;
    }
    verify_iupac(seq_str);

    if (flagTrimMasked()) {
      const auto len = seq_str.length();
      size_t trim_start = 0, trim_end = seq_str.length();
      while (trim_start <= len && bool(islower(seq_str[trim_start]))) {
        trim_start++;
      }
      while (trim_end > 0 && bool(islower(seq_str[trim_end - 1]))) {
        trim_end--;
      }
      seq_str.erase(trim_end);
      seq_str.erase(0, trim_start);
      if (!qual_str.empty()) {
        qual_str.erase(trim_end);
        qual_str.erase(0, trim_start);
      }
    }
    if (flagFoldCase()) {
      std::transform(
        seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
    }
    return true;
  }
  return false;
}

inline std::string
SeqReader::name()
{
  return std::move(name_str);
}

inline std::string
SeqReader::comment()
{
  return std::move(comment_str);
}

inline std::string
SeqReader::seq()
{
  return std::move(seq_str);
}

inline std::string
SeqReader::qual()
{
  return std::move(qual_str);
}

} // namespace btllib

#endif
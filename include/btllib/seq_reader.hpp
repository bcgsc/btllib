#ifndef BTLLIB_SEQ_READER_HPP
#define BTLLIB_SEQ_READER_HPP

#include "data_saveload.hpp"
#include "seq.hpp"
#include "status.hpp"
#include "num_queue.hpp"

#include <atomic>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <string>
#include <thread>

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
  std::atomic<bool> end = false;
  unsigned flags;
  Format format = UNKNOWN; // Format of the input file

  static const size_t RESERVE_SIZE_FOR_STRINGS = 1024;
  static const size_t DETERMINE_FORMAT_CHARS = 2048;
  static const size_t BUFFER_SIZE = 65536;
  static const size_t RECORD_BLOCK_SIZE = 128;
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

  void start_worker();

  bool getline(std::string& line);
  char getsection(std::string& section);

  int read_stage = 0;
  bool read_fasta();
  bool read_fastq();
  bool read_sam();
  bool read_gfa2();

  bool (SeqReader::*read_impl)() = nullptr;

  struct Record {
    std::string name;
    std::string comment;
    std::string seq;
    std::string qual;
  };

  struct RecordBlock {
    size_t num = 0;
    Record records[RECORD_BLOCK_SIZE];
    size_t current = 0;
    size_t count = 0;
  };

  std::thread* worker_thread = nullptr;
  std::string tmp;
  RecordBlock worker_records, ready_records;
  Record *current_worker_record, *current_ready_record;
  InputNumQueue<RecordBlock> queue;
};

inline SeqReader::SeqReader(const std::string& source, int flags)
  : source(source)
  , input(data_load(source))
  , flags(flags)
{
  determine_format();
  tmp.reserve(RESERVE_SIZE_FOR_STRINGS);
  start_worker();
}

inline void
SeqReader::close()
{
  if (!closed) {
    input.close();
    end = true;
    queue.close();
    worker_thread->join();
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
SeqReader::load_buffer()
{
  buffer_start = 0;
  char last = buffer_end > 0 ? buffer[buffer_end - 1] : char(0);
  while ((buffer_end = fread(buffer, 1, BUFFER_SIZE, input)) == 0 &&
         !bool(std::feof(input))) {
  }
  if (bool(std::feof(input)) && !eof_newline_inserted) {
    if (buffer_end < BUFFER_SIZE) {
      if ((buffer_end == 0 && last != '\n') ||
          (buffer_end > 0 && buffer[buffer_end - 1] != '\n')) {
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
  while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
         current < buffer_end) {
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
  while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
         current < buffer_end) {
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

  while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
         current < buffer_end && buffer[current] == '@') {
    while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
           current < buffer_end && buffer[current] != '\n') {
      current++;
    }
    current++;
  }

  int column = 1;
  unsigned char c;
  while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
         current < buffer_end) {
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
  while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
         current < buffer_end) {
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
  check_warning(buffer_end - buffer_start == 1,
                std::string(source) + " is empty.");

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
  char c = char(0);
  for (; buffer_start < buffer_end && (c = buffer[buffer_start]) != '\n';
       ++buffer_start) {
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
  char c = char(0);
  for (; buffer_start < buffer_end && (c = buffer[buffer_start]) != '\n' &&
         c != ' ';
       ++buffer_start) {
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
      if (!getline(tmp)) {
        return false;
      }
      auto pos = tmp.find(' ');
      current_worker_record->name = tmp.substr(1, pos - 1);
      while (pos < tmp.size() && tmp[pos] == ' ') {
        pos++;
      }
      if (pos < tmp.size()) {
        current_worker_record->comment = tmp.substr(pos);
      } else {
        current_worker_record->comment = "";
      }
      ++read_stage;
      tmp.clear();
    }
    case 1: {
      if (!getline(tmp)) {
        return false;
      }
      current_worker_record->seq = std::move(tmp);
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
      if (!getline(tmp)) {
        return false;
      }
      auto pos = tmp.find(' ');
      current_worker_record->name = tmp.substr(1, pos - 1);
      while (pos < tmp.size() && tmp[pos] == ' ') {
        pos++;
      }
      if (pos < tmp.size()) {
        current_worker_record->comment = tmp.substr(pos);
      } else {
        current_worker_record->comment = "";
      }
      ++read_stage;
      tmp.clear();
    }
    case 1: {
      if (!getline(tmp)) {
        return false;
      }
      current_worker_record->seq = std::move(tmp);
      ++read_stage;
      tmp.clear();
    }
    case 2: {
      if (!getline(tmp)) {
        return false;
      }
      ++read_stage;
      tmp.clear();
    }
    case 3: {
      if (!getline(tmp)) {
        return false;
      }
      current_worker_record->qual = std::move(tmp);
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
    if (!getline(tmp)) {
      return false;
    }
    if (tmp.length() > 0 && tmp[0] != '@') {
      size_t pos = 0, pos2 = 0, pos3 = 0;
      pos2 = tmp.find('\t');
      current_worker_record->name = tmp.substr(0, pos2);
      for (int i = 0; i < int(SEQ) - 1; i++) {
        pos = tmp.find('\t', pos + 1);
      }
      pos2 = tmp.find('\t', pos + 1);
      pos3 = tmp.find('\t', pos2 + 1);
      if (pos3 == std::string::npos) {
        pos3 = tmp.length();
      }

      current_worker_record->seq = tmp.substr(pos + 1, pos2 - pos - 1);
      current_worker_record->qual = tmp.substr(pos2 + 1, pos3 - pos2 - 1);

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
    if (!getline(tmp)) {
      return false;
    }
    if (tmp.length() > 0 && tmp[0] == 'S') {
      size_t pos = 0, pos2 = 0;
      pos2 = tmp.find('\t', 1);
      current_worker_record->name = tmp.substr(1, pos2 - 1);
      for (int i = 0; i < int(SEQ) - 1; i++) {
        pos = tmp.find('\t', pos + 1);
      }
      pos2 = tmp.find('\t', pos + 1);
      if (pos2 == std::string::npos) {
        pos2 = tmp.length();
      }

      current_worker_record->seq = tmp.substr(pos + 1, pos2 - pos - 1);

      tmp.clear();
      return true;
    }
    tmp.clear();
    if (buffer_start >= buffer_end) {
      return false;
    }
  }
}

inline void
SeqReader::start_worker() {
  worker_thread = new std::thread([=] () {
    size_t counter = 0;
    for (; buffer_start < buffer_end && !end;) {
      current_worker_record = &(worker_records.records[worker_records.current]);
      while (!(this->*read_impl)()) {
        if (!load_buffer()) {
          break;
        }
      }
      if (buffer_start >= buffer_end) {
        load_buffer();
      }

      if (current_worker_record->seq.empty()) {
        break;
      }
      verify_iupac(current_worker_record->seq);

      if (flagTrimMasked()) {
        const auto len = current_worker_record->seq.length();
        size_t trim_start = 0, trim_end = current_worker_record->seq.length();
        while (trim_start <= len && bool(islower(current_worker_record->seq[trim_start]))) {
          trim_start++;
        }
        while (trim_end > 0 && bool(islower(current_worker_record->seq[trim_end - 1]))) {
          trim_end--;
        }
        current_worker_record->seq.erase(trim_end);
        current_worker_record->seq.erase(0, trim_start);
        if (!current_worker_record->qual.empty()) {
          current_worker_record->qual.erase(trim_end);
          current_worker_record->qual.erase(0, trim_start);
        }
      }
      if (flagFoldCase()) {
        std::transform(
          current_worker_record->seq.begin(), current_worker_record->seq.end(), current_worker_record->seq.begin(), ::toupper);
      }
      worker_records.current++;
      worker_records.count++;
      if (worker_records.current == RECORD_BLOCK_SIZE) {
        worker_records.current = 0;
        worker_records.num = counter++;
        queue.write(worker_records);
        worker_records = RecordBlock();
      }
    }
    end = true;
    worker_records.current = 0;
    worker_records.num = counter++;
    size_t last_count = worker_records.count;
    queue.write(worker_records);
    if (last_count > 0) {
      RecordBlock dummy;
      dummy.num = counter++;
      dummy.count = 0;
      queue.write(dummy);
    }
  });
}

inline bool
SeqReader::read()
{
  if (ready_records.count <= ready_records.current + 1) {
    queue.read(ready_records);
    current_ready_record = &(ready_records.records[0]);
  } else {
    current_ready_record = &(ready_records.records[++ready_records.current]);
  }
  return ready_records.count > 0;
}

inline std::string
SeqReader::name()
{
  return std::move(current_ready_record->name);
}

inline std::string
SeqReader::comment()
{
  return std::move(current_ready_record->comment);
}

inline std::string
SeqReader::seq()
{
  return std::move(current_ready_record->seq);
}

inline std::string
SeqReader::qual()
{
  return std::move(current_ready_record->qual);
}

} // namespace btllib

#endif
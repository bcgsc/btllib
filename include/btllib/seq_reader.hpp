#ifndef BTLLIB_SEQ_READER_HPP
#define BTLLIB_SEQ_READER_HPP

#include "data_saveload.hpp"
#include "num_queue.hpp"
#include "seq.hpp"
#include "status.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <mutex>
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

  SeqReader(const std::string& source_path, int flags = 0);
  ~SeqReader();

  void close();

  bool flagFoldCase() const { return bool(~flags & NO_FOLD_CASE); }
  bool flagTrimMasked() const { return bool(flags & TRIM_MASKED); }

  enum Format
  {
    UNDETERMINED,
    FASTA,
    FASTQ,
    SAM,
    GFA2,
    INVALID
  };

  Format get_format() const { return format; }

  struct Record
  {
    std::string name;
    std::string comment;
    std::string seq;
    std::string qual;

    operator bool() { return !seq.empty(); }
  };

  /** Read operator. */
  Record read();

private:
  const std::string& source_path;
  DataSource source;
  unsigned flags = 0;
  Format format = UNDETERMINED; // Format of the source file
  bool closed = false;
  std::mutex read_mutex;

  static const size_t RESERVE_SIZE_FOR_STRINGS = 1024;
  static const size_t DETERMINE_FORMAT_CHARS = 2048;
  static const size_t BUFFER_SIZE = DETERMINE_FORMAT_CHARS;

  char* buffer = nullptr;
  size_t buffer_start = 0;
  size_t buffer_end = 0;
  bool eof_newline_inserted = false;

  char* line_buffer;
  size_t line_buffer_size = 65536; // NOLINT

  static const size_t RECORD_QUEUE_SIZE = 128;
  static const size_t RECORD_BLOCK_SIZE = 64;

  struct RecordBlock // NOLINT
  {
    size_t num = 0;
    Record records[RECORD_BLOCK_SIZE];
    size_t current = 0;
    size_t count = 0;
  };

  std::thread* worker_thread = nullptr;
  std::atomic<bool> worker_end;
  std::string tmp;
  RecordBlock worker_records, ready_records;
  Record *current_worker_record = nullptr, *current_ready_record = nullptr;
  InputNumQueue<RecordBlock, RECORD_QUEUE_SIZE> queue;

  void determine_format();
  void start_worker();

  bool load_buffer();

  bool is_fasta_buffer();
  bool is_fastq_buffer();
  bool is_sam_buffer();
  bool is_gfa2_buffer();

  bool readline_buffer(std::string& line);
  void readline_file(std::string& line);

  int read_stage = 0;

  bool read_fasta_buffer();
  bool read_fastq_buffer();
  bool read_sam_buffer();
  bool read_gfa2_buffer();

  void read_fasta_file();
  void read_fastq_file();
  void read_sam_file();
  void read_gfa2_file();

  bool (SeqReader::*read_format_buffer_impl)() = nullptr;
  void (SeqReader::*read_format_file_impl)() = nullptr;

  void postprocess();
};

inline SeqReader::SeqReader(const std::string& source_path, int flags)
  : source_path(source_path)
  , source(source_path)
  , flags(flags)
  , worker_end(false)
{
  buffer = new char[BUFFER_SIZE];
  line_buffer = (char*)std::malloc(line_buffer_size); // NOLINT
  tmp.reserve(RESERVE_SIZE_FOR_STRINGS);
  determine_format();
  start_worker();
}

inline SeqReader::~SeqReader()
{
  free(line_buffer); // NOLINT
  close();
}

inline void
SeqReader::close()
{
  if (!closed) {
    source.close();
    worker_end = true;
    queue.close();
    worker_thread->join();
    closed = true;
  }
}

inline bool
SeqReader::load_buffer()
{
  buffer_start = 0;
  char last = buffer_end > 0 ? buffer[buffer_end - 1] : char(0);
  while ((buffer_end = fread(buffer, 1, BUFFER_SIZE, source)) == 0 &&
         !bool(std::feof(source))) {
  }
  if (bool(std::feof(source)) && !eof_newline_inserted) {
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
SeqReader::is_fasta_buffer()
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
SeqReader::is_fastq_buffer()
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
SeqReader::is_sam_buffer()
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
SeqReader::is_gfa2_buffer()
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
                std::string(source_path) + " is empty.");

  if (is_fasta_buffer()) {
    format = Format::FASTA;
    read_format_buffer_impl = &SeqReader::read_fasta_buffer;
    read_format_file_impl = &SeqReader::read_fasta_file;
  } else if (is_fastq_buffer()) {
    format = Format::FASTQ;
    read_format_buffer_impl = &SeqReader::read_fastq_buffer;
    read_format_file_impl = &SeqReader::read_fastq_file;
  } else if (is_sam_buffer()) {
    format = Format::SAM;
    read_format_buffer_impl = &SeqReader::read_sam_buffer;
    read_format_file_impl = &SeqReader::read_sam_file;
  } else if (is_gfa2_buffer()) {
    format = Format::GFA2;
    read_format_buffer_impl = &SeqReader::read_gfa2_buffer;
    read_format_file_impl = &SeqReader::read_gfa2_file;
  } else {
    format = Format::INVALID;
    log_error(std::string(source_path) + " source file is in invalid format!");
    std::exit(EXIT_FAILURE);
  }
}

inline bool
SeqReader::readline_buffer(std::string& line)
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

inline void
SeqReader::readline_file(std::string& line)
{
  size_t p = 0;
  for (;;) {
    line_buffer[line_buffer_size - 1] = char(255); // NOLINT
    fgets(line_buffer + p, int(line_buffer_size - p), source);
    if (line_buffer[line_buffer_size - 1] == 0) {
      p = line_buffer_size - 1;
      line_buffer_size *= 2;
      line_buffer =
        (char*)std::realloc((char*)line_buffer, line_buffer_size); // NOLINT
    } else {
      break;
    }
  }
  line += line_buffer;
  if (line.back() == '\n') {
    line.pop_back();
  }
}

inline bool
SeqReader::read_fasta_buffer()
{
  switch (read_stage) {
    case 0: {
      if (!readline_buffer(tmp)) {
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
      if (!readline_buffer(tmp)) {
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
SeqReader::read_fastq_buffer()
{
  switch (read_stage) {
    case 0: {
      if (!readline_buffer(tmp)) {
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
      if (!readline_buffer(tmp)) {
        return false;
      }
      current_worker_record->seq = std::move(tmp);
      ++read_stage;
      tmp.clear();
    }
    case 2: {
      if (!readline_buffer(tmp)) {
        return false;
      }
      ++read_stage;
      tmp.clear();
    }
    case 3: {
      if (!readline_buffer(tmp)) {
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
SeqReader::read_sam_buffer()
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
    if (!readline_buffer(tmp)) {
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
SeqReader::read_gfa2_buffer()
{
  enum Column
  {
    S = 1,
    ID,
    LEN,
    SEQ
  };
  for (;;) {
    if (!readline_buffer(tmp)) {
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
SeqReader::read_fasta_file()
{
  switch (read_stage) {
    case 0: {
      readline_file(tmp);
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
      readline_file(tmp);
      current_worker_record->seq = std::move(tmp);
      read_stage = 0;
      tmp.clear();
    }
  }
}

inline void
SeqReader::read_fastq_file()
{
  switch (read_stage) {
    case 0: {
      readline_file(tmp);
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
      readline_file(tmp);
      current_worker_record->seq = std::move(tmp);
      ++read_stage;
      tmp.clear();
    }
    case 2: {
      readline_file(tmp);
      ++read_stage;
      tmp.clear();
    }
    case 3: {
      readline_file(tmp);
      current_worker_record->qual = std::move(tmp);
      read_stage = 0;
      tmp.clear();
    }
  }
}

inline void
SeqReader::read_sam_file()
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
    readline_file(tmp);
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
    }
    tmp.clear();
    if (bool(feof(source))) {
      break;
    }
  }
}

inline void
SeqReader::read_gfa2_file()
{
  enum Column
  {
    S = 1,
    ID,
    LEN,
    SEQ
  };
  for (;;) {
    readline_file(tmp);
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
    }
    tmp.clear();
    if (bool(feof(source))) {
      break;
    }
  }
}

inline void
SeqReader::postprocess()
{
  auto& seq = current_worker_record->seq;
  if (flagTrimMasked()) {
    const auto len = seq.length();
    size_t trim_start = 0, trim_end = seq.length();
    while (trim_start <= len && bool(islower(seq[trim_start]))) {
      trim_start++;
    }
    while (trim_end > 0 && bool(islower(seq[trim_end - 1]))) {
      trim_end--;
    }
    seq.erase(trim_end);
    seq.erase(0, trim_start);
    if (!current_worker_record->qual.empty()) {
      current_worker_record->qual.erase(trim_end);
      current_worker_record->qual.erase(0, trim_start);
    }
  }
  if (flagFoldCase()) {
    for (auto& c : seq) {
      c = CAPITALS[unsigned(c)];
      if (!bool(c)) {
        log_error(std::string("A sequence contains invalid IUPAC character: ") +
                  c);
        std::exit(EXIT_FAILURE);
      }
    }
  }
}

inline void
SeqReader::start_worker()
{
  worker_thread = new std::thread([this]() {
    size_t counter = 0;

    // Read from buffer
    for (; buffer_start < buffer_end && !worker_end;) {
      current_worker_record = &(worker_records.records[worker_records.current]);
      if (!(this->*read_format_buffer_impl)() ||
          current_worker_record->seq.empty()) {
        break;
      }
      postprocess();
      worker_records.current++;
      worker_records.count++;
      if (worker_records.current == RECORD_BLOCK_SIZE) {
        worker_records.current = 0;
        worker_records.num = counter++;
        queue.write(worker_records);
        worker_records = RecordBlock();
      }
    }
    delete[] buffer;

    // Read from file (more efficient, less copying involved)
    for (; std::ferror(source) == 0 && std::feof(source) == 0 && !worker_end;) {
      int p = std::fgetc(source);
      if (p == EOF) {
        break;
      }
      std::ungetc(p, source);
      current_worker_record = &(worker_records.records[worker_records.current]);
      (this->*read_format_file_impl)();
      if (current_worker_record->seq.empty()) {
        break;
      }
      postprocess();
      worker_records.current++;
      worker_records.count++;
      if (worker_records.current == RECORD_BLOCK_SIZE) {
        worker_records.current = 0;
        worker_records.num = counter++;
        queue.write(worker_records);
        worker_records = RecordBlock();
      }
    }

    worker_end = true;
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

inline SeqReader::Record
SeqReader::read()
{
  std::unique_lock<std::mutex> lock(read_mutex);
  if (ready_records.count <= ready_records.current + 1) {
    queue.read(ready_records);
    current_ready_record = &(ready_records.records[0]);
  } else {
    current_ready_record = &(ready_records.records[++ready_records.current]);
  }
  if (ready_records.count == 0) {
    close();
  }
  return std::move(*current_ready_record);
}

} // namespace btllib

#endif
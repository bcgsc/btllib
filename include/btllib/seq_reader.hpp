#ifndef BTLLIB_SEQ_READER_HPP
#define BTLLIB_SEQ_READER_HPP

#include "data_saveload.hpp"
#include "index_queue.hpp"
#include "seq.hpp"
#include "status.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cctype>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <stack>
#include <string>
#include <thread>

namespace btllib {

/** Read a FASTA, FASTQ, SAM, or GFA2 file. Threadsafe. */
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

  static const size_t DETERMINE_FORMAT_CHARS = 2048;
  static const size_t BUFFER_SIZE = DETERMINE_FORMAT_CHARS;

  char* buffer = nullptr;
  size_t buffer_start = 0;
  size_t buffer_end = 0;
  bool eof_newline_inserted = false;

  static const size_t RECORD_QUEUE_SIZE = 128;
  static const size_t RECORD_BLOCK_SIZE = 64;

  static const size_t NAME_DEFAULT_CAPACITY = 4096;
  static const size_t COMMENT_DEFAULT_CAPACITY = NAME_DEFAULT_CAPACITY;
  static const size_t SEQ_DEFAULT_CAPACITY = NAME_DEFAULT_CAPACITY;
  static const size_t QUAL_DEFAULT_CAPACITY = SEQ_DEFAULT_CAPACITY;

  static const size_t MAX_SIMULTANEOUS_SEQREADERS = 128;

  struct RecordCString
  {

    RecordCString()
      : name((char*)std::malloc(NAME_DEFAULT_CAPACITY)) // NOLINT
      , name_cap(NAME_DEFAULT_CAPACITY)
      , comment((char*)std::malloc(COMMENT_DEFAULT_CAPACITY)) // NOLINT
      , comment_cap(COMMENT_DEFAULT_CAPACITY)
      , seq((char*)std::malloc(SEQ_DEFAULT_CAPACITY)) // NOLINT
      , seq_cap(SEQ_DEFAULT_CAPACITY)
      , qual((char*)std::malloc(QUAL_DEFAULT_CAPACITY)) // NOLINT
      , qual_cap(QUAL_DEFAULT_CAPACITY)
    {
      name[0] = '\0';
      comment[0] = '\0';
      seq[0] = '\0';
      qual[0] = '\0';
    }

    RecordCString(const RecordCString&) = delete;

    RecordCString(RecordCString&& record) noexcept
    {
      std::swap(name, record.name);
      std::swap(name_cap, record.name_cap);
      std::swap(comment, record.comment);
      std::swap(comment_cap, record.comment_cap);
      std::swap(seq, record.seq);
      std::swap(seq_cap, record.seq_cap);
      std::swap(qual, record.qual);
      std::swap(qual_cap, record.qual_cap);
    }

    RecordCString& operator=(const RecordCString&) = delete;

    RecordCString& operator=(RecordCString&& record) noexcept
    {
      std::swap(name, record.name);
      std::swap(name_cap, record.name_cap);
      std::swap(comment, record.comment);
      std::swap(comment_cap, record.comment_cap);
      std::swap(seq, record.seq);
      std::swap(seq_cap, record.seq_cap);
      std::swap(qual, record.qual);
      std::swap(qual_cap, record.qual_cap);
      return *this;
    }

    ~RecordCString()
    {
      free(name);    // NOLINT
      free(comment); // NOLINT
      free(seq);     // NOLINT
      free(qual);    // NOLINT
    }

    char* name = nullptr;
    size_t name_cap = 0;
    char* comment = nullptr;
    size_t comment_cap = 0;
    char* seq = nullptr;
    size_t seq_cap = 0;
    char* qual = nullptr;
    size_t qual_cap = 0;
  };

  char* tmp = nullptr;
  size_t tmp_cap = 0;

  std::thread* reader_thread = nullptr;
  std::thread* postprocessor_thread = nullptr;
  std::mutex format_mutex;
  std::condition_variable format_cv;
  std::atomic<bool> reader_end;
  RecordCString* reader_record = nullptr;
  IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>
    reader_queue;
  IndexQueueSPMC<Record, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>
    postprocessor_queue;

  thread_local inline static IndexQueueSPMC<Record,
                                            RECORD_QUEUE_SIZE,
                                            RECORD_BLOCK_SIZE>::Block
    ready_records_array[MAX_SIMULTANEOUS_SEQREADERS];
  thread_local inline static Record*
    ready_record_array[MAX_SIMULTANEOUS_SEQREADERS];

  inline static std::stack<unsigned> recycled_ids;
  inline static std::mutex recycled_ids_mutex;
  inline static unsigned last_id = 0;
  unsigned id = 0;

  void determine_format();
  void start_reader();
  void start_postprocessor();

  bool load_buffer();

  bool is_fasta_buffer();
  bool is_fastq_buffer();
  bool is_sam_buffer();
  bool is_gfa2_buffer();

  bool readline_buffer_append(char*& ptr, size_t& cap);
  void readline_file(char*& ptr, size_t& cap);
  void readline_file_append(char*& ptr, size_t& cap);

  int read_stage = 0;

public:
  struct read_fasta_buffer;
  struct read_fastq_buffer;
  struct read_sam_buffer;
  struct read_gfa2_buffer;

  struct read_fasta_transition;
  struct read_fastq_transition;
  struct read_sam_transition;
  struct read_gfa2_transition;

  struct read_fasta_file;
  struct read_fastq_file;
  struct read_sam_file;
  struct read_gfa2_file;

private:
  template<typename F>
  void read_from_buffer(
    F f,
    IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block&
      records,
    size_t& counter);

  template<typename F>
  void read_transition(
    F f,
    IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block&
      records,
    size_t& counter);

  template<typename F>
  void read_from_file(
    F f,
    IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block&
      records,
    size_t& counter);

  void postprocess();
};

inline SeqReader::SeqReader(const std::string& source_path, int flags)
  : source_path(source_path)
  , source(source_path)
  , flags(flags)
  , reader_end(false)
{
  buffer = new char[BUFFER_SIZE];
  tmp = (char*)std::malloc(SEQ_DEFAULT_CAPACITY); // NOLINT
  tmp_cap = SEQ_DEFAULT_CAPACITY;
  tmp[0] = '\0';
  start_postprocessor();
  {
    std::unique_lock<std::mutex> lock(recycled_ids_mutex);
    if (recycled_ids.empty()) {
      id = ++last_id;
    } else {
      id = recycled_ids.top();
      recycled_ids.pop();
    }
  }
  {
    std::unique_lock<std::mutex> lock(format_mutex);
    start_reader();
    format_cv.wait(lock);
  }
}

inline SeqReader::~SeqReader()
{
  {
    std::unique_lock<std::mutex> lock(recycled_ids_mutex);
    recycled_ids.push(id);
  }
  close();
  delete[] buffer;
  free(tmp); // NOLINT
}

inline void
SeqReader::close()
{
  if (!closed) {
    closed = true;
    reader_end = true;
    reader_queue.close();
    postprocessor_queue.close();
    reader_thread->join();
    postprocessor_thread->join();
    source.close();
  }
}

inline bool
SeqReader::load_buffer()
{
  buffer_start = 0;
  char last = buffer_end > 0 ? buffer[buffer_end - 1] : char(0);
  buffer_end = 0;
  do {
    buffer_end +=
      fread(buffer + buffer_end, 1, BUFFER_SIZE - buffer_end, source);
  } while (buffer_end < BUFFER_SIZE && !bool(std::feof(source)));

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
  bool empty = buffer_end - buffer_start == 1;
  check_warning(empty, std::string(source_path) + " is empty.");

  if (empty) {
    return;
  }

  if (is_fasta_buffer()) {
    format = Format::FASTA;
  } else if (is_fastq_buffer()) {
    format = Format::FASTQ;
  } else if (is_sam_buffer()) {
    format = Format::SAM;
  } else if (is_gfa2_buffer()) {
    format = Format::GFA2;
  } else {
    format = Format::INVALID;
    log_error(std::string(source_path) + " source file is in invalid format!");
    std::exit(EXIT_FAILURE);
  }
}

inline bool
SeqReader::readline_buffer_append(char*& ptr, size_t& cap)
{
  char c = char(0);
  size_t size = std::strlen(ptr);
  for (; buffer_start < buffer_end && (c = buffer[buffer_start]) != '\n';
       ++buffer_start) {
    if (size >= cap) {
      cap *= 2;
      ptr = (char*)std::realloc((char*)ptr, cap); // NOLINT
    }
    ptr[size++] = c;
  }
  if (size >= cap) {
    cap *= 2;
    ptr = (char*)std::realloc((char*)ptr, cap); // NOLINT
  }
  ptr[size++] = '\0';
  if (c == '\n') {
    ++buffer_start;
    return true;
  }
  return false;
}

inline void
SeqReader::readline_file(char*& ptr, size_t& cap)
{
  getline(&ptr, &cap, source);
}

inline void
SeqReader::readline_file_append(char*& ptr, size_t& cap)
{
  readline_file(tmp, tmp_cap);
  size_t tmp_size = std::strlen(tmp);
  size_t ptr_size = std::strlen(ptr);
  while (tmp_size + ptr_size > cap) {
    cap *= 2;
    ptr = (char*)std::realloc((char*)ptr, cap); // NOLINT
  }
  strcat(ptr, tmp); // NOLINT
}

// NOLINTNEXTLINE
#define READ_FASTA_NAME_COMMENT(READLINE_SECTION, END_SECTION)                 \
  READLINE_SECTION                                                             \
  size_t tmp_size = strlen(seq_reader.tmp);                                    \
  while (tmp_size > seq_reader.reader_record->name_cap) {                      \
    seq_reader.reader_record->name_cap *= 2;                                   \
    seq_reader.reader_record->name =                                           \
      (char*)std::realloc((char*)seq_reader.reader_record->name,               \
                          seq_reader.reader_record->name_cap);                 \
  }                                                                            \
  while (tmp_size > seq_reader.reader_record->comment_cap) {                   \
    seq_reader.reader_record->comment_cap *= 2;                                \
    seq_reader.reader_record->comment = /* NOLINTNEXTLINE */                   \
      (char*)std::realloc((char*)seq_reader.reader_record->comment,            \
                          seq_reader.reader_record->comment_cap);              \
  }                                                                            \
  bool in_name = true;                                                         \
  char c;                                                                      \
  int i, j;                                                                    \
  for (i = 0, j = 0, c = seq_reader.tmp[1]; c != '\0';                         \
       c = seq_reader.tmp[++i + 1]) {                                          \
    if (in_name) {                                                             \
      if (c == ' ') {                                                          \
        seq_reader.reader_record->name[i] = '\0';                              \
        in_name = false;                                                       \
      } else {                                                                 \
        seq_reader.reader_record->name[i] = c;                                 \
      }                                                                        \
    } else {                                                                   \
      seq_reader.reader_record->comment[j++] = c;                              \
    }                                                                          \
  }                                                                            \
  if (in_name) {                                                               \
    seq_reader.reader_record->name[++i] = '\0';                                \
  }                                                                            \
  seq_reader.reader_record->comment[j] = '\0';                                 \
  END_SECTION

// NOLINTNEXTLINE
#define READ_SAM(READLINE_SECTION, MIDEND_SECTION, END_SECTION)                \
  enum Column                                                                  \
  {                                                                            \
    QNAME = 1,                                                                 \
    FLAG,                                                                      \
    RNAME,                                                                     \
    POS,                                                                       \
    MAPQ,                                                                      \
    CIGAR,                                                                     \
    RNEXT,                                                                     \
    PNEXT,                                                                     \
    TLEN,                                                                      \
    SEQ,                                                                       \
    QUAL                                                                       \
  };                                                                           \
  for (;;) {                                                                   \
    READLINE_SECTION                                                           \
    std::string tmp_string = seq_reader.tmp;                                   \
    if (tmp_string.length() > 0 && tmp_string[0] != '@') {                     \
      size_t pos = 0, pos2 = 0, pos3 = 0;                                      \
      pos2 = tmp_string.find('\t');                                            \
      while (tmp_string.size() > seq_reader.reader_record->name_cap) {         \
        seq_reader.reader_record->name_cap *= 2;                               \
        seq_reader.reader_record->name =                                       \
          (char*)std::realloc((char*)seq_reader.reader_record->name,           \
                              seq_reader.reader_record->name_cap);             \
      }                                                                        \
      strcpy(seq_reader.reader_record->name,                                   \
             tmp_string.substr(0, pos2).c_str());                              \
      for (int i = 0; i < int(SEQ) - 1; i++) {                                 \
        pos = tmp_string.find('\t', pos + 1);                                  \
      }                                                                        \
      pos2 = tmp_string.find('\t', pos + 1);                                   \
      pos3 = tmp_string.find('\t', pos2 + 1);                                  \
      if (pos3 == std::string::npos) {                                         \
        pos3 = tmp_string.length();                                            \
      }                                                                        \
      while (tmp_string.size() > seq_reader.reader_record->seq_cap) {          \
        seq_reader.reader_record->seq_cap *= 2;                                \
        seq_reader.reader_record->seq =                                        \
          (char*)std::realloc((char*)seq_reader.reader_record->seq,            \
                              seq_reader.reader_record->seq_cap);              \
      }                                                                        \
      while (tmp_string.size() > seq_reader.reader_record->qual_cap) {         \
        seq_reader.reader_record->qual_cap *= 2;                               \
        seq_reader.reader_record->qual =                                       \
          (char*)std::realloc((char*)seq_reader.reader_record->qual,           \
                              seq_reader.reader_record->qual_cap);             \
      }                                                                        \
      strcpy(seq_reader.reader_record->seq,                                    \
             tmp_string.substr(pos + 1, pos2 - pos - 1).c_str());              \
      strcpy(seq_reader.reader_record->qual,                                   \
             tmp_string.substr(pos2 + 1, pos3 - pos2 - 1).c_str());            \
      MIDEND_SECTION                                                           \
    }                                                                          \
    seq_reader.tmp[0] = '\0';                                                  \
    END_SECTION                                                                \
  }

// NOLINTNEXTLINE
#define READ_GFA2(READLINE_SECTION, MIDEND_SECTION, END_SECTION)               \
  enum Column                                                                  \
  {                                                                            \
    S = 1,                                                                     \
    ID,                                                                        \
    LEN,                                                                       \
    SEQ                                                                        \
  };                                                                           \
  for (;;) {                                                                   \
    READLINE_SECTION                                                           \
    std::string tmp_string = seq_reader.tmp;                                   \
    if (tmp_string.length() > 0 && tmp_string[0] == 'S') {                     \
      size_t pos = 0, pos2 = 0;                                                \
      pos2 = tmp_string.find('\t', 1);                                         \
      while (tmp_string.size() > seq_reader.reader_record->name_cap) {         \
        seq_reader.reader_record->name_cap *= 2;                               \
        seq_reader.reader_record->name =                                       \
          (char*)std::realloc((char*)seq_reader.reader_record->name,           \
                              seq_reader.reader_record->name_cap);             \
      }                                                                        \
      strcpy(seq_reader.reader_record->name,                                   \
             tmp_string.substr(1, pos2 - 1).c_str());                          \
      for (int i = 0; i < int(SEQ) - 1; i++) {                                 \
        pos = tmp_string.find('\t', pos + 1);                                  \
      }                                                                        \
      pos2 = tmp_string.find('\t', pos + 1);                                   \
      if (pos2 == std::string::npos) {                                         \
        pos2 = tmp_string.length();                                            \
      }                                                                        \
      while (tmp_string.size() > seq_reader.reader_record->seq_cap) {          \
        seq_reader.reader_record->seq_cap *= 2;                                \
        seq_reader.reader_record->seq =                                        \
          (char*)std::realloc((char*)seq_reader.reader_record->seq,            \
                              seq_reader.reader_record->seq_cap);              \
      }                                                                        \
      strcpy(seq_reader.reader_record->seq,                                    \
             tmp_string.substr(pos + 1, pos2 - pos - 1).c_str());              \
      MIDEND_SECTION                                                           \
    }                                                                          \
    seq_reader.tmp[0] = '\0';                                                  \
    END_SECTION                                                                \
  }

struct SeqReader::read_fasta_buffer
{
  bool operator()(SeqReader& seq_reader)
  {
    switch (seq_reader.read_stage) {
      case 0: {
        READ_FASTA_NAME_COMMENT(                        // NOLINT
          if (!seq_reader.readline_buffer_append(       // NOLINT
                seq_reader.tmp,                         // NOLINT
                seq_reader.tmp_cap)) { return false; }, // NOLINT
          ++seq_reader.read_stage;                      // NOLINT
          seq_reader.tmp[0] = '\0';)                    // NOLINT
      }
      case 1: {
        if (!seq_reader.readline_buffer_append(
              seq_reader.reader_record->seq,
              seq_reader.reader_record->seq_cap)) {
          return false;
        }
        seq_reader.read_stage = 0;
        return true;
      }
    }
    return false;
  }
};

struct SeqReader::read_fastq_buffer
{
  bool operator()(SeqReader& seq_reader)
  {
    switch (seq_reader.read_stage) {
      case 0: {
        READ_FASTA_NAME_COMMENT(                        // NOLINT
          if (!seq_reader.readline_buffer_append(       // NOLINT
                seq_reader.tmp,                         // NOLINT
                seq_reader.tmp_cap)) { return false; }, // NOLINT
          ++seq_reader.read_stage;                      // NOLINT
          seq_reader.tmp[0] = '\0';)                    // NOLINT
      }
      case 1: {
        if (!seq_reader.readline_buffer_append(
              seq_reader.reader_record->seq,
              seq_reader.reader_record->seq_cap)) {
          return false;
        }
        ++seq_reader.read_stage;
      }
      case 2: {
        if (!seq_reader.readline_buffer_append(seq_reader.tmp,
                                               seq_reader.tmp_cap)) {
          return false;
        }
        ++seq_reader.read_stage;
        seq_reader.tmp[0] = '\0';
      }
      case 3: {
        if (!seq_reader.readline_buffer_append(
              seq_reader.reader_record->qual,
              seq_reader.reader_record->qual_cap)) {
          return false;
        }
        seq_reader.read_stage = 0;
        return true;
      }
    }
    return false;
  }
};

struct SeqReader::read_sam_buffer
{
  bool operator()(SeqReader& seq_reader)
  {
    READ_SAM(                                       // NOLINT
      if (!seq_reader.readline_buffer_append(       // NOLINT
            seq_reader.tmp,                         // NOLINT
            seq_reader.tmp_cap)) { return false; }, // NOLINT
      seq_reader.tmp[0] = '\0';                     // NOLINT
      return true;                                  // NOLINT
      ,
      if (seq_reader.buffer_start >= seq_reader.buffer_end) {
        return false;
      }) // NOLINT
  }
};

struct SeqReader::read_gfa2_buffer
{
  bool operator()(SeqReader& seq_reader)
  {
    READ_GFA2(                                // NOLINT
      if (!seq_reader.readline_buffer_append( // NOLINT
            seq_reader.tmp,
            seq_reader.tmp_cap)) { return false; }, // NOLINT
      seq_reader.tmp[0] = '\0';                     // NOLINT
      return true;                                  // NOLINT
      ,
      if (seq_reader.buffer_start >= seq_reader.buffer_end) {
        return false;
      }) // NOLINT
  }
};

struct SeqReader::read_fasta_transition
{
  void operator()(SeqReader& seq_reader)
  {
    switch (seq_reader.read_stage) {
      case 0: {
        READ_FASTA_NAME_COMMENT(                               // NOLINT
          seq_reader.readline_file_append(seq_reader.tmp,      // NOLINT
                                          seq_reader.tmp_cap); // NOLINT
          , ++seq_reader.read_stage;                           // NOLINT
          seq_reader.tmp[0] = '\0';)                           // NOLINT
      }
      case 1: {
        seq_reader.readline_file_append(seq_reader.reader_record->name,
                                        seq_reader.reader_record->name_cap);
        seq_reader.read_stage = 0;
      }
    }
  }
};

struct SeqReader::read_fastq_transition
{
  void operator()(SeqReader& seq_reader)
  {
    switch (seq_reader.read_stage) {
      case 0: {
        READ_FASTA_NAME_COMMENT(                               // NOLINT
          seq_reader.readline_file_append(seq_reader.tmp,      // NOLINT
                                          seq_reader.tmp_cap); // NOLINT
          , ++seq_reader.read_stage;                           // NOLINT
          seq_reader.tmp[0] = '\0';)                           // NOLINT
      }
      case 1: {
        seq_reader.readline_file_append(seq_reader.reader_record->seq,
                                        seq_reader.reader_record->seq_cap);
        ++seq_reader.read_stage;
      }
      case 2: {
        seq_reader.readline_file_append(seq_reader.tmp, seq_reader.tmp_cap);
        ++seq_reader.read_stage;
        seq_reader.tmp[0] = '\0';
      }
      case 3: {
        seq_reader.readline_file_append(seq_reader.reader_record->qual,
                                        seq_reader.reader_record->qual_cap);
        seq_reader.read_stage = 0;
      }
    }
  }
};

struct SeqReader::read_sam_transition
{
  void operator()(SeqReader& seq_reader)
  {
    READ_SAM(                                              // NOLINT
      seq_reader.readline_file_append(seq_reader.tmp,      // NOLINT
                                      seq_reader.tmp_cap); // NOLINT
      , , if (bool(feof(seq_reader.source))) { break; })   // NOLINT
  }
};

struct SeqReader::read_gfa2_transition
{
  void operator()(SeqReader& seq_reader)
  {
    READ_GFA2(                                             // NOLINT
      seq_reader.readline_file_append(seq_reader.tmp,      // NOLINT
                                      seq_reader.tmp_cap); // NOLINT
      , , if (bool(feof(seq_reader.source))) { break; })   // NOLINT
  }
};

struct SeqReader::read_fasta_file
{
  void operator()(SeqReader& seq_reader)
  {
    READ_FASTA_NAME_COMMENT(                                           // NOLINT
      seq_reader.readline_file(seq_reader.tmp, seq_reader.tmp_cap);, ) // NOLINT
    seq_reader.readline_file(seq_reader.reader_record->seq,
                             seq_reader.reader_record->seq_cap);
  }
};

struct SeqReader::read_fastq_file
{
  void operator()(SeqReader& seq_reader)
  {
    READ_FASTA_NAME_COMMENT(                                           // NOLINT
      seq_reader.readline_file(seq_reader.tmp, seq_reader.tmp_cap);, ) // NOLINT
    seq_reader.readline_file(seq_reader.reader_record->seq,
                             seq_reader.reader_record->seq_cap);
    seq_reader.readline_file(seq_reader.tmp, seq_reader.tmp_cap);
    seq_reader.readline_file(seq_reader.reader_record->qual,
                             seq_reader.reader_record->qual_cap);
  }
};

struct SeqReader::read_sam_file
{
  void operator()(SeqReader& seq_reader)
  {
    READ_SAM(                                                       // NOLINT
      seq_reader.readline_file(seq_reader.tmp, seq_reader.tmp_cap); // NOLINT
      , , if (bool(feof(seq_reader.source))) { break; })            // NOLINT
  }
};

struct SeqReader::read_gfa2_file
{
  void operator()(SeqReader& seq_reader)
  {
    READ_GFA2(                                                      // NOLINT
      seq_reader.readline_file(seq_reader.tmp, seq_reader.tmp_cap); // NOLINT
      , , if (bool(feof(seq_reader.source))) { break; })            // NOLINT
  }
};

template<typename F>
inline void
SeqReader::read_from_buffer(
  F f,
  IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block&
    records,
  size_t& counter)
{
  for (; buffer_start < buffer_end && !reader_end;) {
    reader_record = &(records.data[records.current]);
    if (!f(*this) || reader_record->seq[0] == '\0') {
      break;
    }
    records.current++;
    records.count++;
    if (records.current == RECORD_BLOCK_SIZE) {
      records.current = 0;
      records.index = counter++;
      reader_queue.write(records);
      records = IndexQueueSPSC<RecordCString,
                               RECORD_QUEUE_SIZE,
                               RECORD_BLOCK_SIZE>::Block();
    }
  }
}

template<typename F>
inline void
SeqReader::read_transition(
  F f,
  IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block&
    records,
  size_t& counter)
{
  if (std::ferror(source) == 0 && std::feof(source) == 0 && !reader_end) {
    int p = std::fgetc(source);
    if (p != EOF) {
      std::ungetc(p, source);
      reader_record = &(records.data[records.current]);
      f(*this);
      if (reader_record->seq[0] != '\0') {
        records.current++;
        records.count++;
        if (records.current == RECORD_BLOCK_SIZE) {
          records.current = 0;
          records.index = counter++;
          reader_queue.write(records);
          records = IndexQueueSPSC<RecordCString,
                                   RECORD_QUEUE_SIZE,
                                   RECORD_BLOCK_SIZE>::Block();
        }
      }
    }
  }
}

template<typename F>
inline void
SeqReader::read_from_file(
  F f,
  IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block&
    records,
  size_t& counter)
{
  for (; std::ferror(source) == 0 && std::feof(source) == 0 && !reader_end;) {
    int p = std::fgetc(source);
    if (p == EOF) {
      break;
    }
    std::ungetc(p, source);
    reader_record = &(records.data[records.current]);
    f(*this);
    if (reader_record->seq[0] == '\0') {
      break;
    }
    records.current++;
    records.count++;
    if (records.current == RECORD_BLOCK_SIZE) {
      records.current = 0;
      records.index = counter++;
      reader_queue.write(records);
      records = IndexQueueSPSC<RecordCString,
                               RECORD_QUEUE_SIZE,
                               RECORD_BLOCK_SIZE>::Block();
    }
  }
}

inline void
SeqReader::start_reader()
{
  reader_thread = new std::thread([this]() {
    {
      std::unique_lock<std::mutex> lock(format_mutex);
      determine_format();
      format_cv.notify_all();
    }

    size_t counter = 0;
    IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block
      records;
    switch (format) {
      case FASTA: {
        read_from_buffer(read_fasta_buffer(), records, counter);
        read_transition(read_fasta_transition(), records, counter);
        read_from_file(read_fasta_file(), records, counter);
        break;
      }
      case FASTQ: {
        read_from_buffer(read_fastq_buffer(), records, counter);
        read_transition(read_fastq_transition(), records, counter);
        read_from_file(read_fastq_file(), records, counter);
        break;
      }
      case SAM: {
        read_from_buffer(read_sam_buffer(), records, counter);
        read_transition(read_sam_transition(), records, counter);
        read_from_file(read_sam_file(), records, counter);
        break;
      }
      case GFA2: {
        read_from_buffer(read_gfa2_buffer(), records, counter);
        read_transition(read_gfa2_transition(), records, counter);
        read_from_file(read_gfa2_file(), records, counter);
        break;
      }
      default: {
        break;
      }
    }

    reader_end = true;
    records.current = 0;
    records.index = counter++;
    size_t last_count = records.count;
    reader_queue.write(records);
    if (last_count > 0) {
      IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block
        dummy;
      dummy.index = counter++;
      dummy.current = 0;
      dummy.count = 0;
      reader_queue.write(dummy);
    }
  });
}

inline void
SeqReader::start_postprocessor()
{
  postprocessor_thread = new std::thread([this]() {
    IndexQueueSPSC<RecordCString, RECORD_QUEUE_SIZE, RECORD_BLOCK_SIZE>::Block
      records_reader;
    decltype(postprocessor_queue)::Block records_postprocessor;
    for (;;) {
      reader_queue.read(records_reader);
      for (size_t i = 0; i < records_reader.count; i++) {
        records_postprocessor.data[i].name = records_reader.data[i].name;
        records_postprocessor.data[i].comment = records_reader.data[i].comment;
        records_postprocessor.data[i].seq = records_reader.data[i].seq;
        records_postprocessor.data[i].qual = records_reader.data[i].qual;
        auto& name = records_postprocessor.data[i].name;
        auto& comment = records_postprocessor.data[i].comment;
        auto& seq = records_postprocessor.data[i].seq;
        auto& qual = records_postprocessor.data[i].qual;
        if (!name.empty() && name.back() == '\n') {
          name.pop_back();
        }
        if (!comment.empty() && comment.back() == '\n') {
          comment.pop_back();
        }
        if (!seq.empty() && seq.back() == '\n') {
          seq.pop_back();
        }
        if (!qual.empty() && qual.back() == '\n') {
          qual.pop_back();
        }
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
          if (!qual.empty()) {
            qual.erase(trim_end);
            qual.erase(0, trim_start);
          }
        }
        if (flagFoldCase()) {
          for (auto& c : seq) {
            char old = c;
            c = CAPITALS[unsigned(c)];
            if (!bool(c)) {
              log_error(std::string("A sequence contains invalid "
                                    "IUPAC character: ") +
                        old);
              std::exit(EXIT_FAILURE);
            }
          }
        }
        records_reader.data[i].name[0] = '\0';
        records_reader.data[i].comment[0] = '\0';
        records_reader.data[i].seq[0] = '\0';
        records_reader.data[i].qual[0] = '\0';
      }
      records_postprocessor.count = records_reader.count;
      records_postprocessor.current = records_reader.current;
      records_postprocessor.index = records_reader.index;
      if (records_postprocessor.count == 0) {
        postprocessor_queue.write(records_postprocessor);
        break;
      }
      postprocessor_queue.write(records_postprocessor);
    }
  });
}

inline SeqReader::Record
SeqReader::read()
{
  auto& ready_records = ready_records_array[id];
  auto& ready_record = ready_record_array[id];
  if (ready_records.count <= ready_records.current) {
    postprocessor_queue.read(ready_records);
    if (ready_records.count <= ready_records.current) {
      close();
      return Record();
    }
  }
  ready_record = &(ready_records.data[ready_records.current++]);
  return std::move(*ready_record);
}

} // namespace btllib

#endif
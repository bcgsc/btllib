#ifndef BTLLIB_SEQ_READER_SAM_MODULE_HPP
#define BTLLIB_SEQ_READER_SAM_MODULE_HPP

#include "cstring.hpp"
#include "seq.hpp"

namespace btllib {

/// @cond HIDDEN_SYMBOLS
class SeqReaderSamModule
{

private:
  friend class SeqReader;

  enum class Stage
  {
    HEADER,
    ALIGNMENTS
  };

  Stage stage = Stage::HEADER;
  CString tmp;

  bool buffer_valid(const char* buffer, const size_t size) const;
  template<typename ReaderType, typename RecordType>
  bool read_buffer(ReaderType& reader, RecordType& record);
  template<typename ReaderType, typename RecordType>
  void read_transition(ReaderType& reader, RecordType& record);
  template<typename ReaderType, typename RecordType>
  void read_file(ReaderType& reader, RecordType& record);
};

inline bool
SeqReaderSamModule::buffer_valid(const char* buffer, const size_t size) const
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

  while (current < size && buffer[current] == '@') {
    while (current < size && buffer[current] != '\n') {
      current++;
    }
    current++;
  }

  Column column = QNAME;
  unsigned char c;
  while (current < size) {
    c = buffer[current];
    if (c == '\n') {
      break;
    }
    if (c == '\t') {
      if (current > 0 && !bool(std::isspace(buffer[current - 1]))) {
        column = Column(int(column) + 1);
      } else {
        return false;
      }
    } else {
      switch (column) {
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

  return current >= size || column >= QUAL;
}

// NOLINTNEXTLINE
#define READ_SAM(READLINE_SECTION, STOP_SECTION, RETURN_SECTION)               \
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
                                                                               \
  READLINE_SECTION                                                             \
  if (stage == Stage::HEADER) {                                                \
    for (;;) {                                                                 \
      if (tmp[0] != '@') {                                                     \
        stage = Stage::ALIGNMENTS;                                             \
        break;                                                                 \
      }                                                                        \
      STOP_SECTION                                                             \
      tmp.clear();                                                             \
      READLINE_SECTION                                                         \
    }                                                                          \
  }                                                                            \
                                                                               \
  std::string tmp_string(tmp);                                                 \
  tmp.clear();                                                                 \
                                                                               \
  const size_t qname_end = tmp_string.find('\t');                              \
  size_t pos = qname_end;                                                      \
  for (int i = 1; i < int(SEQ) - 1; i++) {                                     \
    pos = tmp_string.find('\t', pos + 1);                                      \
  }                                                                            \
  const size_t seq_start = pos + 1;                                            \
  const size_t seq_end = tmp_string.find('\t', seq_start + 1);                 \
  const size_t qual_start = seq_end + 1;                                       \
  size_t qual_end = tmp_string.find('\t', qual_start + 1);                     \
  if (qual_end == std::string::npos) {                                         \
    qual_end = tmp_string.length();                                            \
  }                                                                            \
                                                                               \
  if (qname_end + 1 > record.header.s_cap) {                                   \
    record.header.change_cap(qname_end + 1);                                   \
  }                                                                            \
  if (seq_end - seq_start + 1 > record.seq.s_cap) {                            \
    record.seq.change_cap(seq_end - seq_start + 1);                            \
  }                                                                            \
  if (qual_end - qual_start + 1 > record.qual.s_cap) {                         \
    record.qual.change_cap(qual_end - qual_start + 1);                         \
  }                                                                            \
                                                                               \
  record.header = tmp_string.substr(0, qname_end);                             \
  record.seq = tmp_string.substr(seq_start, seq_end - seq_start);              \
  record.qual = tmp_string.substr(qual_start, qual_end - qual_start);          \
  RETURN_SECTION

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderSamModule::read_buffer(ReaderType& reader, RecordType& record)
{
  READ_SAM(                                                          // NOLINT
    if (!reader.readline_buffer_append(tmp)) { return false; },      // NOLINT
    if (reader.buffer.start >= reader.buffer.end) { return false; }, // NOLINT
    return true;)                                                    // NOLINT
}

template<typename ReaderType, typename RecordType>
inline void
SeqReaderSamModule::read_transition(ReaderType& reader, RecordType& record)
{
  READ_SAM(                                       // NOLINT
    reader.readline_file_append(tmp);             // NOLINT
    ,                                             // NOLINT
    if (bool(feof(reader.source))) { return; }, ) // NOLINT
}

template<typename ReaderType, typename RecordType>
inline void
SeqReaderSamModule::read_file(ReaderType& reader, RecordType& record)
{
  READ_SAM(                                       // NOLINT
    reader.readline_file(tmp);                    // NOLINT
    ,                                             // NOLINT
    if (bool(feof(reader.source))) { return; }, ) // NOLINT
}
/// @endcond

} // namespace btllib

#endif
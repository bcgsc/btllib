#ifndef BTLLIB_SEQ_READER_GFA2_MODULE_HPP
#define BTLLIB_SEQ_READER_GFA2_MODULE_HPP

#include "cstring.hpp"
#include "seq.hpp"

namespace btllib {

/// @cond HIDDEN_SYMBOLS
class SeqReaderGfa2Module
{

private:
  friend class SeqReader;

  enum class Stage
  {
    HEADER,
    SEQ,
    SEP,
    QUAL
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

// // NOLINTNEXTLINE
// #define READ_GFA2(READLINE_SECTION, MIDEND_SECTION, END_SECTION) \
//   enum Column \
//   { \
//     S = 1, \
//     ID, \
//     LEN, \
//     SEQ \
//   }; \
//   for (;;) { \
//     READLINE_SECTION \
//     std::string tmp_string = seq_reader.tmp.s; \
//     if (tmp_string.length() > 0 && tmp_string[0] == 'S') { \
//       size_t pos = 0, pos2 = 0; \
//       pos2 = tmp_string.find('\t', 1); \
//       if (tmp_string.size() + 1 > seq_reader.reader_record->header.s_cap) { \
//         seq_reader.reader_record->header.s_cap = tmp_string.size() + 1; \
//         seq_reader.reader_record->header.s = \
//           (char*)std::realloc((char*)(seq_reader.reader_record->header.s), \
//                               seq_reader.reader_record->header.s_cap); \
//       } \
//       seq_reader.reader_record->header = tmp_string.substr(1, pos2 - 1); \
//       for (int i = 0; i < int(SEQ) - 1; i++) { \
//         pos = tmp_string.find('\t', pos + 1); \
//       } \
//       pos2 = tmp_string.find('\t', pos + 1); \
//       if (pos2 == std::string::npos) { \
//         pos2 = tmp_string.length(); \
//       } \
//       if (tmp_string.size() + 1 > seq_reader.reader_record->seq.s_cap) { \
//         seq_reader.reader_record->seq.s_cap = tmp_string.size() + 1; \
//         seq_reader.reader_record->seq.s = \
//           (char*)std::realloc((char*)(seq_reader.reader_record->seq.s), \
//                               seq_reader.reader_record->seq.s_cap); \
//       } \
//       seq_reader.reader_record->seq = \
//         tmp_string.substr(pos + 1, pos2 - pos - 1); \
//       MIDEND_SECTION \
//     } \
//     seq_reader.tmp.clear(); \
//     END_SECTION \
//   }

inline bool
SeqReaderGfa2Module::buffer_valid(const char* buffer, const size_t size) const
{
//   const unsigned char specs[] = { 'H', 'S', 'F', 'E', 'G', 'O', 'U' };

//   enum State
//   {
//     IN_ID,
//     IN_ID_TAB,
//     IN_REST,
//     IN_IGNORED
//   };

//   auto is_a_spec = [&](unsigned char c) {
//     bool found = false;
//     for (unsigned char spec : specs) {
//       if (c == spec) {
//         found = true;
//         break;
//       }
//     }
//     return found;
//   };

//   State state = is_a_spec(buffer[0]) ? IN_ID : IN_IGNORED;
//   bool has_id = false;
//   size_t current = buffer_start;
//   unsigned char c;
//   while (current < buffer_start + DETERMINE_FORMAT_CHARS &&
//          current < buffer_end) {
//     c = buffer[current];
//     switch (state) {
//       case IN_ID:
//         if (!is_a_spec(c)) {
//           return false;
//         }
//         has_id = true;
//         state = IN_ID_TAB;
//         break;
//       case IN_ID_TAB:
//         if (c != '\t') {
//           return false;
//         }
//         state = IN_REST;
//         break;
//       case IN_REST:
//         if (c == '\n') {
//           if (current + 1 < buffer_end) {
//             state = is_a_spec(buffer[current + 1]) ? IN_ID : IN_IGNORED;
//           }
//         }
//         break;
//       case IN_IGNORED:
//         if (c == '\n') {
//           if (current + 1 < buffer_end) {
//             state = is_a_spec(buffer[current + 1]) ? IN_ID : IN_IGNORED;
//           }
//         }
//         break;
//       default:
//         break;
//     }
//     current++;
//   }

//   return has_id;
  (void)buffer;
  (void)size;
  return false;
}

template<typename ReaderType, typename RecordType>
inline bool
SeqReaderGfa2Module::read_buffer(ReaderType& reader, RecordType& record)
{
//     READ_GFA2(                                  // NOLINT
//       if (!seq_reader.readline_buffer_append(   // NOLINT
//             seq_reader.tmp)) { return false; }, // NOLINT
//       seq_reader.tmp.clear();                   // NOLINT
//       return true;                              // NOLINT
//       ,
//       if (seq_reader.buffer_start >= seq_reader.buffer_end) {
//         return false;
//       }) // NOLINT
  (void)reader;
  (void)record;
  return false;
}

template<typename ReaderType, typename RecordType>
inline void
SeqReaderGfa2Module::read_transition(ReaderType& reader, RecordType& record)
{
//     READ_GFA2(                                           // NOLINT
//       seq_reader.readline_file_append(seq_reader.tmp);   // NOLINT
//       , , if (bool(feof(seq_reader.source))) { break; }) // NOLINT
  (void)reader;
  (void)record;
}

template<typename ReaderType, typename RecordType>
inline void
SeqReaderGfa2Module::read_file(ReaderType& reader, RecordType& record)
{
//     READ_GFA2(                                           // NOLINT
//       seq_reader.readline_file(seq_reader.tmp);          // NOLINT
//       , , if (bool(feof(seq_reader.source))) { break; }) // NOLINT
  (void)reader;
  (void)record;
}
/// @endcond

} // namespace btllib

#endif
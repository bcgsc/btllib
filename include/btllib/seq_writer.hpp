#ifndef BTLLIB_SEQ_WRITER_HPP
#define BTLLIB_SEQ_WRITER_HPP

#include "data_saveload.hpp"
#include "seq.hpp"

#include <cstdio>
#include <string>

namespace btllib {

static inline ssize_t
fd_write(int fd, const void* buf, size_t count)
{
  return write(fd, buf, count);
}

class SeqWriter
{

public:
  enum Format
  {
    FASTA,
    FASTQ
  };

  SeqWriter(const std::string& sink_path,
            Format format = FASTA,
            bool append = false);

  void close();

  void write(const std::string& name,
             const std::string& comment,
             const std::string& seq,
             const std::string& qual);

private:
  const std::string sink_path;
  DataSink sink;
  bool closed;
  Format format;
  char headerchar;
};

inline SeqWriter::SeqWriter(const std::string& sink_path,
                            Format format,
                            bool append)
  : sink_path(sink_path)
  , sink(sink_path, append)
  , closed(false)
  , format(format)
  , headerchar(format == FASTA ? '>' : '@')
{}

inline void
SeqWriter::close()
{
  if (!closed) {
    sink.close();
    closed = true;
  }
}

inline void
SeqWriter::write(const std::string& name,
                 const std::string& comment,
                 const std::string& seq,
                 const std::string& qual)
{
  check_error(seq.empty(), "Attempted to write empty sequence.");
  for (const auto& c : seq) {
    if (!bool(COMPLEMENTS[unsigned(c)])) {
      log_error(std::string("A sequence contains invalid IUPAC character: ") +
                c);
      std::exit(EXIT_FAILURE);
    }
  }

  std::string output;
  output.reserve(1 + name.size() + 1 + comment.size() + 1 + seq.size() + 3 +
                 qual.size() + 1);
  output += headerchar;
  if (!name.empty()) {
    output += name;
  }
  if (!comment.empty()) {
    output += " ";
    output += comment;
    output += '\n';
  }

  output += seq;
  output += '\n';

  if (format == FASTQ) {
    check_error(seq.size() != qual.size(),
                "Quality must be the same length as sequence.");
    output += "+\n";
    output += qual;
    output += '\n';
  }

  fd_write(sink, output.c_str(), output.size());
}

} // namespace btllib

#endif
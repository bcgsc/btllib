#ifndef BTLLIB_SEQ_WRITER_HPP
#define BTLLIB_SEQ_WRITER_HPP

#include "data_saveload.hpp"
#include "seq.hpp"

#include <cstdio>
#include <string>

namespace btllib {

class SeqWriter
{

public:
  enum Format
  {
    FASTA,
    FASTQ
  };

  SeqWriter(const std::string& sink_path, Format format, bool append = false);

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

  fwrite(&headerchar, 1, 1, sink);
  if (!name.empty()) {
    fwrite(name.c_str(), 1, name.size(), sink);
  }
  if (!comment.empty()) {
    fwrite(" ", 1, 1, sink);
    fwrite(comment.c_str(), 1, comment.size(), sink);
    fwrite("\n", 1, 1, sink);
  }
  fwrite(seq.c_str(), 1, seq.size(), sink);
  fwrite("\n", 1, 1, sink);
  if (format == FASTQ) {
    check_error(seq.size() != qual.size(),
                "Quality must be the same length as sequence.");
    fwrite("+\n", 1, 2, sink);
    fwrite(qual.c_str(), 1, qual.size(), sink);
    fwrite("\n", 1, 1, sink);
  }
}

} // namespace btllib

#endif
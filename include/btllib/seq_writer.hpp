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

  SeqWriter(const std::string& sink, Format format, bool append = false);

  void close();

  void write(const std::string& name,
             const std::string& comment,
             const std::string& seq,
             const std::string& qual);

private:
  const std::string sink;
  DataSink output;
  bool closed;
  Format format;
  char headerchar;
};

inline SeqWriter::SeqWriter(const std::string& sink, Format format, bool append)
  : sink(sink)
  , output(data_save(sink, append))
  , closed(false)
  , format(format)
  , headerchar(format == FASTA ? '>' : '@')
{}

inline void
SeqWriter::close()
{
  if (!closed) {
    output.close();
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
  verify_iupac(seq);

  fwrite(&headerchar, 1, 1, output);
  if (!name.empty()) {
    fwrite(name.c_str(), 1, name.size(), output);
  }
  if (!comment.empty()) {
    fwrite(" ", 1, 1, output);
    fwrite(comment.c_str(), 1, comment.size(), output);
    fwrite("\n", 1, 1, output);
  }
  fwrite(seq.c_str(), 1, seq.size(), output);
  fwrite("\n", 1, 1, output);
  if (format == FASTQ) {
    check_error(seq.size() != qual.size(),
                "Quality must be the same length as sequence.");
    fwrite("+\n", 1, 2, output);
    fwrite(qual.c_str(), 1, qual.size(), output);
    fwrite("\n", 1, 1, output);
  }
}

} // namespace btllib

#endif
#ifndef BTL_SEQUENCE_READER_HPP
#define BTL_SEQUENCE_READER_HPP

#include "status.hpp"
#include "sequence.hpp"

#include <cassert>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <iostream>

namespace btl {

/** Read a FASTA, FASTQ, or SAM file. */
class SequenceReader
{
public:
    enum Flags
    {
        /** Fold lower-case characters to upper-case. */
        FOLD_CASE = 0,
        NO_FOLD_CASE = 1,
        /** Convert to standard quality. */
        NO_CONVERT_QUALITY = 0,
        CONVERT_QUALITY = 2,
        /** Trim masked (lower case) characters from the ends of
	        * sequences. */
        NO_TRIM_MASKED = 0,
        TRIM_MASKED = 4
    };

    SequenceReader(const char* input_path, int flags = 0);

  	bool flagFoldCase() const { return ~flags & NO_FOLD_CASE; }
  	bool flagConvertQual() const { return flags & CONVERT_QUALITY; }
    bool flagTrimMasked() const { return flags & TRIM_MASKED; }

    enum Format {
      FASTA,
      FASTQ,
      SAM
    };

    Format get_format() const { return format; }

    /** Return whether this stream is at end-of-file. */
    bool eof() const { return is.eof(); };

    /** Return true if failbit or badbit of stream is set. */
    bool fail() const { return is.fail(); };

    /** Return whether this stream is good. */
    operator const void*() const { return is ? this : NULL; }

    /** Return the next character of this stream. */
    int peek() { return is.peek(); }

    /** Interface for manipulators. */
    SequenceReader& operator<<(std::istream& (*f)(std::istream&));

    /** Read operator. */
    SequenceReader& operator>>(std::string& seq);

    /** Quality of the last read sequence. */
    std::string get_qual();

private:
    const char* input_path;
    std::ifstream ifs; // If reading from a file
    std::istream& is; // Ref to the input stream (file or stdin)
    int flags;
    Format format; // Format of the input file

    void determine_format();
};

inline SequenceReader::SequenceReader(const char* input_path, int flags)
    : input_path(input_path)
    , ifs(input_path)
    , is(strcmp(input_path, "-") == 0 ? std::cin : ifs)
    , flags(flags)
{
    if (strcmp(input_path, "-") != 0) {
      check_stream(ifs, input_path);
    }
    check_warning(is.eof(), input_path, " is empty.");
    determine_format();
}

inline void SequenceReader::determine_format() {
    while (is.peek() == '@') {
      is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    if (is.peek() == '>') {
      format = FASTA;
    } else {
      while
    }
    is.seekg(0); 
}

inline SequenceReader& SequenceReader::operator<<(std::istream& (*f)(std::istream&))
{
  	f(is);
  	return *this;
}

inline SequenceReader& SequenceReader::operator>>(std::string& seq)
{
	  std::string id, comment, qual;
	  char anchor;
	  return *this;
}

inline std::string get_qual() {

}

} // namespace btl

#endif
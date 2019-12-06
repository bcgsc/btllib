#ifndef BTL_SEQUENCE_READER_HPP
#define BTL_SEQUENCE_READER_HPP

#include "status.hpp"
#include "sequence.hpp"

#include <cassert>
#include <cstring>
#include <cctype>
#include <fstream>
#include <sstream>
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
        UNKNOWN,
        FASTA,
        FASTQ,
        SAM,
        INVALID
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

    static const std::streamsize DETERMINE_FORMAT_MAX_READ_CHARS = 4096;

    void determine_format();
};

inline SequenceReader::SequenceReader(const char* input_path, int flags)
    : input_path(input_path)
    , ifs(input_path)
    , is(strcmp(input_path, "-") == 0 ? std::cin : ifs)
    , flags(flags)
    , format(Format::UNKNOWN)
{
    if (strcmp(input_path, "-") != 0) {
      check_stream(ifs, input_path);
    }
    check_warning(is.eof(), input_path, " is empty.");
    determine_format();
}

inline bool is_fasta(const char* input, size_t n) {
    int current = 0;

}

inline bool is_fastq(const char* input, size_t n) {

}

inline bool is_sam(const char* input, size_t n) {
    int current = 0;

    while (current < n && input[current] == '@') { current++; }

    int column = 1;
    char c;
    while (current < n) {
        c = input[current];
        if (c == '\n') {
            break;
        } else if (c == '\t') {
            if (current > 0 && !std::isspace(input[current - 1])) {
                column++;
            } else {
                return false;
            }
        } else {
            switch (column) {
                case 1: if (std::isspace(c)) { return false; } break;
                case 2: if (!std::isdigit(c)) { return false; } break;
                case 3: if (std::isspace(c)) { return false; } break;
                case 4: if (!std::isdigit(c)) { return false; } break;
                case 5: if (!std::isdigit(c)) { return false; } break;
                case 6: if (std::isspace(c)) { return false; } break;
                case 7: if (std::isspace(c)) { return false; } break;
                case 8: if (!std::isdigit(c)) { return false; } break;
                case 9: if (!std::isdigit(c)) { return false; } break;
                case 10: if (!COMPLEMENTS[c]) { return false; } break;
                case 11: if (std::isspace(c)) { return false; } break;
                default: break;
            }
        }
        current++;
    }
    if (current >= n || column >= 11) {
        return true;
    }

    return false;
}

inline void SequenceReader::determine_format() {
    std::ifstream ifs(input_path);
    std::istream& is(strcmp(input_path, "-") == 0 ? std::cin : ifs);

    char buffer[DETERMINE_FORMAT_MAX_READ_CHARS];
    is.read(buffer, DETERMINE_FORMAT_MAX_READ_CHARS);

    if (is_fasta(buffer, DETERMINE_FORMAT_MAX_READ_CHARS)) {
        format = Format::FASTA;
    } else if (is_fastq(buffer, DETERMINE_FORMAT_MAX_READ_CHARS)) {
        format = Format::FASTQ;
    } else if (is_sam(buffer, DETERMINE_FORMAT_MAX_READ_CHARS)) {
        format == Format::SAM;
    } else {
        format == Format::INVALID;
        raise_error(input_path, " input file is in invalid format!");
    }
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
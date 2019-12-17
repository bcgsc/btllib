#ifndef BTL_SEQ_READER_HPP
#define BTL_SEQ_READER_HPP

#include "status.hpp"
#include "seq.hpp"

#include <algorithm>
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

    SeqReader(const char* input_path, int flags = 0);

  	bool flagFoldCase() const { return ~flags & NO_FOLD_CASE; }
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

    /** Read operator. */
    SeqReader& read(std::string& seq);
    SeqReader& operator>>(std::string& seq) { return read(seq); }

    /** Interface for manipulators. */
    SeqReader& manip(std::istream& (*f)(std::istream&)); 
    SeqReader& operator<<(std::istream& (*f)(std::istream&)) { return manip(f); }

    /** Quality of the last read sequence. */
    std::string get_qual();

private:
    const char* input_path;
    std::ifstream ifs; // If reading from a file
    std::istream& is; // Ref to the input stream (file or stdin)
    int flags;
    Format format; // Format of the input file

    void (SeqReader::*read_impl)(std::string& seq);
    std::string qual;
    std::string tmp;

    static const std::streamsize DETERMINE_FORMAT_MAX_READ_CHARS = 4096;

    static bool is_fasta(const char* input, size_t n);    
    static bool is_fastq(const char* input, size_t n);
    static bool is_sam(const char* input, size_t n);

    void determine_format();

    void read_fasta(std::string& seq);
    void read_fastq(std::string& seq);
    void read_sam(std::string& seq);
};

inline SeqReader::SeqReader(const char* input_path, int flags)
    : input_path(input_path)
    , ifs(input_path)
    , is(strcmp(input_path, "-") == 0 ? std::cin : ifs)
    , flags(flags)
    , format(Format::UNKNOWN)
    , read_impl(nullptr)
{
    if (strcmp(input_path, "-") != 0) {
      check_stream(ifs, input_path);
    }
    check_warning(is.eof(), std::string(input_path) + " is empty.");
    determine_format();
}

inline bool SeqReader::is_fasta(const char* input, size_t n) {
    size_t current = 0;
    unsigned char c;
    enum State { IN_HEADER_1, IN_HEADER_2, IN_SEQ };
    State state = IN_HEADER_1;
    while (current < n) {
        c = input[current];
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
                } else if (!COMPLEMENTS[c]) {
                    return false;
                }
                break;
        }
        current++;
    }
    return true;
}

inline bool SeqReader::is_fastq(const char* input, size_t n) {
    size_t current = 0;
    unsigned char c;
    enum State { IN_HEADER_1, IN_HEADER_2, IN_SEQ, IN_PLUS_1, IN_PLUS_2, IN_QUAL };
    State state = IN_HEADER_1;
    while (current < n) {
        c = input[current];
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
                    state = IN_PLUS_1;
                } else if (!COMPLEMENTS[c]) {
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

inline bool SeqReader::is_sam(const char* input, size_t n) {
    size_t current = 0;

    while (current < n && input[current] == '@') {
        while (current < n && input[current] != '\n') {
            current++;
        }
        current++;
    }

    int column = 1;
    unsigned char c;
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

inline void SeqReader::determine_format() {
    std::ifstream ifs(input_path);
    std::istream& is(strcmp(input_path, "-") == 0 ? std::cin : ifs);

    char buffer[DETERMINE_FORMAT_MAX_READ_CHARS];
    is.read(buffer, DETERMINE_FORMAT_MAX_READ_CHARS);

    if (is_fasta(buffer, is.gcount())) {
        format = Format::FASTA;
        read_impl = &SeqReader::read_fasta;
    } else if (is_fastq(buffer, is.gcount())) {
        format = Format::FASTQ;
        read_impl= &SeqReader::read_fastq;
    } else if (is_sam(buffer, is.gcount())) {
        format = Format::SAM;
        read_impl= &SeqReader::read_sam;
    } else {
        format = Format::INVALID;
        log_error(std::string(input_path) + " input file is in invalid format!");
        std::exit(EXIT_FAILURE);
    }
}

void SeqReader::read_fasta(std::string& seq) {
    std::getline(is, tmp);
    std::getline(is, seq);
}

void SeqReader::read_fastq(std::string& seq) {
    std::getline(is, tmp);
    std::getline(is, seq);
    std::getline(is, tmp);
    std::getline(is, qual);
}

void SeqReader::read_sam(std::string& seq) {
    for (;;) {
        std::getline(is, tmp);
        if (tmp.length() > 0 && tmp[0] != '@') {
            size_t pos = 0, pos2 = 0, pos3 = 0;
            for (int i = 0; i < 9; i++) {
                pos = tmp.find('\t', pos + 1);
            }
            pos2 = tmp.find('\t', pos + 1);
            pos3 = tmp.find('\t', pos2 + 1);
            if (pos3 == std::string::npos) { pos3 = tmp.length(); }

            seq = tmp.substr(pos + 1, pos2 - pos - 1);
            qual = tmp.substr(pos2 + 1, pos3 - pos2 - 1);

            break;
        }
        if (is.eof()) {
            break;
        }
    }
}

inline SeqReader& SeqReader::manip(std::istream& (*f)(std::istream&))
{
  	f(is);
  	return *this;
}

inline SeqReader& SeqReader::read(std::string& seq)
{
    (this->*read_impl)(seq);
    if (flagTrimMasked()) {
        const auto len = seq.length();
        size_t trim_start = 0, trim_end = seq.length();
        while (trim_start <= len && islower(seq[trim_start])) {
            trim_start++;
        }
        while (trim_end > 0 && islower(seq[trim_end - 1])) {
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
        std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    }
	return *this;
}

inline std::string SeqReader::get_qual() {
    return qual;
}

} // namespace btl

#endif
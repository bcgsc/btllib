#ifndef BTL_SEQUENCE_READER_HPP
#define BTL_SEQUENCE_READER_HPP

#include "check.hpp"

#include <cassert>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <iostream>

namespace btl {

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

  SequenceReader(const char* filepath, int flags = 0);

	bool flagFoldCase() const { return ~flags & NO_FOLD_CASE; }
	bool flagConvertQual() const { return flags & CONVERT_QUALITY; }

  /** Split the fasta file into nsections and seek to the start
		* of section. */
	void split(unsigned section, unsigned nsections);

	/** Return whether this stream is at end-of-file. */
	bool eof() const { return is.eof(); };

	/** Return true if failbit or badbit of stream is set. */
	bool fail() const { return is.fail(); };

  /** Return whether this stream is good. */
	operator const void*() const { return is ? this : NULL; }

	/** Return the next character of this stream. */
	int peek() { return is.peek(); }

  /** Interface for manipulators. */
	SequenceReader& operator>>(std::istream& (*f)(std::istream&));

	SequenceReader& operator>>(std::string& seq);

private:
  const char* filepath;
  std::ifstream ifs;
  std::istream& is;
  int flags;

  enum Format {
    FASTA,
    FASTQ,
    SAM
  };

  /** Format of the input file. */
  Format format;

  void determine_format();

  static const char COMPLEMENTS[256];
  static const char CAPITALS[256];
};

inline SequenceReader::SequenceReader(const char* filepath, int flags)
  : filepath(filepath)
  , ifs(filepath)
  , is(strcmp(filepath, "-") == 0 ? std::cin : ifs)
  , flags(flags)
{
  if (strcmp(filepath, "-") != 0) {
    check_stream(ifs, filepath);
  }
  check_warning(is.eof(), filepath, " is empty.");
  determine_format();
}

inline void SequenceReader::split(unsigned section, unsigned nsections) {
  assert(nsections >= section);
	assert(section > 0);
	assert(strcmp(filepath, "-") != 0);

	if (nsections == 1) { return; }
  
	is.seekg(0, std::ios_base::end);
	std::streampos length = is.tellg();
	std::streampos start = length * (section - 1) / nsections;
	std::streampos end = length * section / nsections;
	if (end < length) {
		is.seekg(end);
		if (is.peek() == '>')
			end += 1;
	}

	is.seekg(start);
	if (start > 0) {
		is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		is.ignore(std::numeric_limits<std::streamsize>::max(), '>');
    check_warning(is.eof(), filepath, ": ", section, ": no contigs in this section");
		is.putback('>');
	}
  check_stream(is, filepath);
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

inline SequenceReader& SequenceReader::operator>>(std::istream& (*f)(std::istream&))
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

const inline char SequenceReader::COMPLEMENTS[256] = {
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//       !    "    #    $    %    &    '    (    )    *    +    ,    -    .    /
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , '-' , '.',  0 ,

//  0    1    2    3    4    5    6    7    8    9    :    ;    <    =    >    ?
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//  @    A    B    C    D    E    F    G    H    I    J    K    L    M    N    O
    0 , 'T', 'V', 'G', 'H',  0 ,  0 , 'C', 'D',  0 ,  0 , 'M',  0 , 'K', 'N',  0 , 

//  P    Q    R    S    T    U    V    W    X    Y    Z    [    \    ]    ^    _
    0 ,  0 , 'Y', 'S', 'A', 'U', 'B', 'W',  0 , 'R',  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//  `    a    b    c    d    e    f    g    h    i    j    k    l    m    n    o
    0 , 't', 'v', 'g', 'h',  0 ,  0 , 'c', 'd',  0 ,  0 , 'm',  0 , 'k', 'n',  0 ,

//  p    q    r    s    t    u    v    w    x    y    z    {    |    }    ~   DEL
    0 ,  0 , 'y', 's', 'a', 'u', 'b', 'w',  0 , 'r',  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0
};

const inline char SequenceReader::CAPITALS[256] = {
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//       !    "    #    $    %    &    '    (    )    *    +    ,    -    .    /
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , '-' , '.',  0 ,

//  0    1    2    3    4    5    6    7    8    9    :    ;    <    =    >    ?
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,

//  @    A    B    C    D    E    F    G    H    I    J    K    L    M    N    O
    0 , 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 

//  P    Q    R    S    T    U    V    W    X    Y    Z    [    \    ]    ^    _
   'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',  0 ,  0 ,  0 ,  0 ,  0 ,

//  `    a    b    c    d    e    f    g    h    i    j    k    l    m    n    o
    0 , 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',

//  p    q    r    s    t    u    v    w    x    y    z    {    |    }    ~   DEL
   'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',  0 ,  0 ,  0 ,  0 ,  0 ,

    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
    0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0
};

} // namespace btl

#endif
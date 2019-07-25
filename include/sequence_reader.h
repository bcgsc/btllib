#ifndef SEQUENCE_READER_H
#define SEQUENCE_READER_H

#include <iostream>
#include <fstream>
#include <cstring>

namespace btl {

class SequenceReader {

public:

    SequenceReader(const char* filepath);

private:

    std::ifstream ifs;
    std::istream& is;

};

SequenceReader::SequenceReader(const char* filepath):
    ifs(filepath),
    is(strcmp(filepath, "-") == 0 ? std::cin : ifs)
{
    if (strcmp(filepath, "-") != 0) {
    }
}

} // namespace btl

#endif
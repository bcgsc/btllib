// #ifndef BTL_SEQUENCE_READER_H
// #define BTL_SEQUENCE_READER_H

// #include <iostream>
// #include <fstream>
// #include <cstring>
// #include <thread>

// namespace btl {

// class SequenceReader {

// public:

//     SequenceReader(const char* filepath);

// private:

//     static const size_t BUFFER_SIZE = 1024 * 1024;

//     const char buffer[BUFFER_SIZE];

//     std::ifstream ifs;
//     std::istream& is;
//     std::thread worker;

// };

// SequenceReader::SequenceReader(const char* filepath):
//     ifs(filepath),
//     is(strcmp(filepath, "-") == 0 ? std::cin : ifs)
// {
//     if (strcmp(filepath, "-") != 0) {
//     }
// }

// } // namespace btl

// #endif
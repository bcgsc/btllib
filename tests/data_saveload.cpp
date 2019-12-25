#include "../include/btllib/data_saveload.hpp"

#include <fstream>

int main() {
    // Data loading is already tested in `seq_reader.cpp`, so saving primarily is tested here.

    std::string txt = "data_saveload test";
    std::string line;

    // Test .gz
    const char* gz_filename = "test.gz";

    std::ofstream gz_ostream(gz_filename);
    gz_ostream << txt;
    gz_ostream.close();

    std::ifstream gz_istream(gz_filename);
    getline(gz_istream, line);
    gz_istream.close();
    assert(line == txt);

    // Test .xz
    const char* xz_filename = "test.xz";
    
    std::ofstream xz_ostream(xz_filename);
    xz_ostream << txt;
    xz_ostream.close();

    std::ifstream xz_istream(xz_filename);
    getline(xz_istream, line);
    xz_istream.close();
    assert(line == txt);

    return 0;
}
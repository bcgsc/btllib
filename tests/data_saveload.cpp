#include "../include/btllib/data_saveload.hpp"

#include <fstream>
#include <chrono>
#include <thread>
#include <cstring>

int main() {
    // Data loading is already tested in `seq_reader.cpp`, so saving primarily is tested here.

    const char* txt = "data_saveload test";
    char* line = new char[128];
    size_t line_len;

    // Test .gz
    const char* gz_filename = "test.gz";

    auto gz_sink = btllib::DataSink(gz_filename, false);
    fwrite(txt, strlen(txt), 1, gz_sink);
    gz_sink.close();

    auto gz_source = btllib::DataSource(gz_filename);
    getline(&line, &line_len, gz_source);
    gz_source.close();
    assert(strcmp(line, txt) == 0);

    // Test .xz
    const char* xz_filename = "test.xz";

    auto xz_sink = btllib::DataSink(xz_filename, false);
    fwrite(txt, strlen(txt), 1, xz_sink);
    xz_sink.close();

    auto xz_source = btllib::DataSource(xz_filename);
    getline(&line, &line_len, xz_source);
    xz_source.close();
    assert(strcmp(line, txt) == 0);

    return 0;
}
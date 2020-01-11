#include "../include/btllib/data_saveload.hpp"

#include <fstream>
#include <chrono>
#include <thread>
#include <cstring>
#include <cstdio>

int main() {
    // Data loading is already tested in `seq_reader.cpp`, so saving primarily is tested here.

    const char* txt = "data_saveload test";
    char* line = new char[128];
    size_t line_len;

    // Test .gz
    const char* gz_filename = "test.gz";

    auto gz_sink = btllib::DataSink(gz_filename, false);
    auto gz_sink_file = fdopen(gz_sink, "w");
    fwrite(txt, strlen(txt), 1, gz_sink_file);
    fclose(gz_sink_file);
    gz_sink.close();

    auto gz_source = btllib::DataSource(gz_filename);
    auto gz_source_file = fdopen(gz_source, "r");
    getline(&line, &line_len, gz_source_file);
    fclose(gz_sink_file);
    gz_source.close();
    assert(strcmp(line, txt) == 0);

    std::remove("test.gz");

    // Test .xz
    const char* xz_filename = "test.xz";

    auto xz_sink = btllib::DataSink(xz_filename, false);
    auto xz_sink_file = fdopen(xz_sink, "w");
    fwrite(txt, strlen(txt), 1, xz_sink_file);
    fclose(xz_sink_file);
    xz_sink.close();

    auto xz_source = btllib::DataSource(xz_filename);
    auto xz_source_file = fdopen(xz_source, "r");
    getline(&line, &line_len, xz_source_file);
    fclose(xz_source_file);
    xz_source.close();
    assert(strcmp(line, txt) == 0);

    std::remove("test.xz");

    return 0;
}
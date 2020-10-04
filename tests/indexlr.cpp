#include "../include/btllib/indexlr.hpp"

#include <fstream>
#include <sstream>

int
main()
{
  btllib::Indexlr indexlr("../tests/indexlr.fa", 100, 5, 0, 1);
  decltype(indexlr)::Record record;
  std::stringstream ss;
  int i = 0;
  while ((record = indexlr.get_minimizers())) {
    if (i > 0) {
      ss << '\n';
    }
    ss << record.id << '\t';
    int j = 0;
    for (const auto& min : record.minimizers) {
      if (j > 0) {
        ss << ' ';
      }
      ss << min.hash2;
      j++;
    }
    i++;
  }

  std::ifstream correct_output_file("../tests/indexlr.fa.correct");
  std::string correct_output;
  correct_output_file.seekg(0, std::ios::end);
  correct_output.reserve(correct_output_file.tellg());
  correct_output_file.seekg(0, std::ios::beg);

  correct_output.assign(std::istreambuf_iterator<char>(correct_output_file),
                        std::istreambuf_iterator<char>());

  assert(ss.str() == correct_output);

  btllib::Indexlr indexlr2("../tests/indexlr.fq",
                           100,
                           5,
                           btllib::Indexlr::Flag::ID |
                             btllib::Indexlr::Flag::BX |
                             btllib::Indexlr::Flag::SEQ,
                           1);

  ss.str("");
  i = 0;
  while ((record = indexlr2.get_minimizers())) {
    if (i > 0) {
      ss << '\n';
    }
    ss << record.id << '\t' << record.barcode << '\t';
    int j = 0;
    for (const auto& min : record.minimizers) {
      if (j > 0) {
        ss << ' ';
      }
      ss << min.hash2 << ':' << min.pos << ':' << min.seq;
      j++;
    }
    i++;
  }

  std::ifstream correct_output_file2("../tests/indexlr.fq.correct");
  correct_output_file2.seekg(0, std::ios::end);
  correct_output.reserve(correct_output_file2.tellg());
  correct_output_file2.seekg(0, std::ios::beg);

  correct_output.assign(std::istreambuf_iterator<char>(correct_output_file2),
                        std::istreambuf_iterator<char>());

  assert(ss.str() == correct_output);

  return 0;
}
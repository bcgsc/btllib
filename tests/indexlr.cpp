#include "../include/btllib/indexlr.hpp"

#include <fstream>
#include <sstream>

int
main()
{
  btllib::Indexlr indexlr("../tests/indexlr.fa", 100, 5, 0);
  btllib::Indexlr indexlr2("../tests/indexlr.fq",
                           100,
                           5,
                           btllib::Indexlr::Flag::ID |
                             btllib::Indexlr::Flag::BX |
                             btllib::Indexlr::Flag::SEQ);

  std::ifstream correct_output_file("../tests/indexlr.fa.correct");
  std::string correct_output;
  correct_output_file.seekg(0, std::ios::end);
  correct_output.reserve(correct_output_file.tellg());
  correct_output_file.seekg(0, std::ios::beg);

  correct_output.assign(std::istreambuf_iterator<char>(correct_output_file),
                        std::istreambuf_iterator<char>());

  std::ifstream correct_output_file2("../tests/indexlr.fq.correct");
  std::string correct_output2;
  correct_output_file2.seekg(0, std::ios::end);
  correct_output2.reserve(correct_output_file2.tellg());
  correct_output_file2.seekg(0, std::ios::beg);

  correct_output2.assign(std::istreambuf_iterator<char>(correct_output_file2),
                         std::istreambuf_iterator<char>());

  std::stringstream ss;
  std::stringstream ss2;

  decltype(indexlr)::Record record;
  bool success_indexlr = false, success_indexlr2 = false;
  for (int i = 0;; i++) {
    if ((success_indexlr = (record = indexlr.get_minimizers()))) {
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
    }
    if ((success_indexlr2 = (record = indexlr2.get_minimizers()))) {
      if (i > 0) {
        ss2 << '\n';
      }
      ss2 << record.id << '\t' << record.barcode << '\t';
      int j = 0;
      for (const auto& min : record.minimizers) {
        if (j > 0) {
          ss2 << ' ';
        }
        ss2 << min.hash2 << ':' << min.pos << ':' << min.seq;
        j++;
      }
    }
    if (!success_indexlr && !success_indexlr2) {
      break;
    }
  }

  assert(ss.str() == correct_output);
  assert(ss2.str() == correct_output2);

  return 0;
}
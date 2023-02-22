#include "btllib/indexlr.hpp"
#include "btllib/bloom_filter.hpp"
#include "btllib/util.hpp"
#include "helpers.hpp"

#include <fstream>
#include <sstream>

int
main()
{
  std::cerr << "Testing on empty file" << std::endl;
  btllib::Indexlr indexlr_empty(btllib::get_dirname(__FILE__) + "/empty.fa",
                                20,
                                50,
                                btllib::Indexlr::Flag::LONG_MODE);
  int i = 0;
  for (auto minimizers : indexlr_empty) {
    i++;
  }
  TEST_ASSERT_EQ(i, 0);

  std::cerr << "Testing without Bloom filters" << std::endl;
  btllib::Indexlr indexlr(btllib::get_dirname(__FILE__) + "/indexlr.fa",
                          100,
                          5,
                          btllib::Indexlr::Flag::SHORT_MODE);
  btllib::Indexlr indexlr2(btllib::get_dirname(__FILE__) + "/indexlr.fq",
                           100,
                           5,
                           btllib::Indexlr::Flag::BX |
                             btllib::Indexlr::Flag::SEQ |
                             btllib::Indexlr::Flag::SHORT_MODE);
  btllib::Indexlr indexlr3(btllib::get_dirname(__FILE__) + "/indexlr.stlfr.fq",
                           75,
                           5,
                           btllib::Indexlr::Flag::BX |
                             btllib::Indexlr::Flag::SHORT_MODE |
                             btllib::Indexlr::Flag::NO_ID);
  
  std::cerr << "Testing quality filtering" << std::endl;
  btllib::Indexlr indexlr4(btllib::get_dirname(__FILE__) + "/indexlr.quality.fq",
                            50,
                            20,
                            size_t(25),
                             btllib::Indexlr::Flag::SHORT_MODE);
                                                          

  std::ifstream correct_output_file(btllib::get_dirname(__FILE__) +
                                    "/indexlr.fa.correct");
  std::string correct_output;
  correct_output_file.seekg(0, std::ios::end);
  correct_output.reserve(correct_output_file.tellg());
  correct_output_file.seekg(0, std::ios::beg);

  correct_output.assign(std::istreambuf_iterator<char>(correct_output_file),
                        std::istreambuf_iterator<char>());

  std::ifstream correct_output_file2(btllib::get_dirname(__FILE__) +
                                     "/indexlr.fq.correct");
  std::string correct_output2;
  correct_output_file2.seekg(0, std::ios::end);
  correct_output2.reserve(correct_output_file2.tellg());
  correct_output_file2.seekg(0, std::ios::beg);

  correct_output2.assign(std::istreambuf_iterator<char>(correct_output_file2),
                         std::istreambuf_iterator<char>());

  std::ifstream correct_output_file3(btllib::get_dirname(__FILE__) +
                                     "/indexlr.stlfr.fq.correct");
  std::string correct_output3;
  correct_output_file3.seekg(0, std::ios::end);
  correct_output3.reserve(correct_output_file3.tellg());
  correct_output_file3.seekg(0, std::ios::beg);

  correct_output3.assign(std::istreambuf_iterator<char>(correct_output_file3),
                         std::istreambuf_iterator<char>());

  std::ifstream correct_output_file4(btllib::get_dirname(__FILE__) +
                                     "/indexlr.quality.fq.correct");
  std::string correct_output4;
  correct_output_file4.seekg(0, std::ios::end);
  correct_output4.reserve(correct_output_file4.tellg());
  correct_output_file4.seekg(0, std::ios::beg);

  correct_output4.assign(std::istreambuf_iterator<char>(correct_output_file4),
                         std::istreambuf_iterator<char>());

  std::stringstream ss;
  std::stringstream ss2;
  std::stringstream ss3;
  std::stringstream ss4;

  decltype(indexlr)::Record record;
  bool success_indexlr = false, success_indexlr2 = false, success_indexlr3 = false, success_indexlr4 = false;
  for (int i = 0;; i++) {
    if ((success_indexlr = (record = indexlr.read()))) {
      if (i > 0) {
        ss << '\n';
      }
      ss << record.id << '\t';
      int j = 0;
      for (const auto& min : record.minimizers) {
        if (j > 0) {
          ss << ' ';
        }
        ss << min.out_hash;
        j++;
      }
    }
    if ((success_indexlr2 = (record = indexlr2.read()))) {
      if (i > 0) {
        ss2 << '\n';
      }
      ss2 << record.id << '\t' << record.barcode << '\t';
      int j = 0;
      for (const auto& min : record.minimizers) {
        if (j > 0) {
          ss2 << ' ';
        }
        ss2 << min.out_hash << ':' << min.pos << ':'
            << (min.forward ? '+' : '-') << ':' << min.seq;
        j++;
      }
    }
    if ((success_indexlr3 = (record = indexlr3.read()))) {
      if (i > 0) {
        ss3 << '\n';
      }
      ss3 << record.barcode << '\t';
      int j = 0;
      for (const auto& min : record.minimizers) {
        if (j > 0) {
          ss3 << ' ';
        }
        ss3 << min.out_hash;
        j++;
      }
    }
    if ((success_indexlr4 = (record = indexlr4.read()))) {
      if (i > 0) {
        ss4 << '\n';
      }
      ss4 << record.id << '\t';
      int j = 0;
      for (const auto& min : record.minimizers) {
        if (j > 0) {
          ss4 << ' ';
        }
        ss4 << min.out_hash;
        j++;
      }
    }
    if (!success_indexlr && !success_indexlr2 && !success_indexlr3 && !success_indexlr4) {
      break;
    }
  }

  ss << std::endl;
  ss2 << std::endl;
  ss3 << std::endl;
  ss4 << std::endl;

  TEST_ASSERT_EQ(ss.str(), correct_output);
  TEST_ASSERT_EQ(ss2.str(), correct_output2);
  TEST_ASSERT_EQ(ss3.str(), correct_output3);
  TEST_ASSERT_EQ(ss4.str(), correct_output4);


  std::cerr << "Testing with Bloom filters" << std::endl;
  btllib::BloomFilter filter_in_bf(1024 * 1024 * 32, 1);
  btllib::BloomFilter filter_out_bf(1024 * 1024 * 32, 1);

  std::vector<uint64_t> filter_in_hashes = { 1315163655624994337ULL,
                                             4261937130627716230ULL,
                                             6055236352205909654ULL };
  std::vector<uint64_t> filter_out_hashes = { 54854220342742384ULL,
                                              10545593919277017579ULL,
                                              6459583435485319281ULL };

  for (const auto h : filter_in_hashes) {
    filter_in_bf.insert({ h });
  }
  for (const auto h : filter_out_hashes) {
    filter_out_bf.insert({ h });
  }

  btllib::Indexlr indexlr5(btllib::get_dirname(__FILE__) + "/indexlr.fq",
                           100,
                           5,
                           btllib::Indexlr::Flag::FILTER_IN |
                             btllib::Indexlr::Flag::LONG_MODE,
                           3,
                           true,
                           filter_in_bf);
  size_t mins_found = 0;
  while ((record = indexlr5.read())) {
    for (const auto& min : record.minimizers) {
      bool found = false;
      for (const auto h : filter_in_hashes) {
        if (min.min_hash == h) {
          found = true;
          break;
        }
      }
      TEST_ASSERT(found);
      mins_found++;
    }
  }
  TEST_ASSERT_GE(mins_found, filter_in_hashes.size());

  btllib::Indexlr indexlr6(btllib::get_dirname(__FILE__) + "/indexlr.fq",
                           100,
                           5,
                           btllib::Indexlr::Flag::FILTER_OUT |
                             btllib::Indexlr::Flag::LONG_MODE,
                           3,
                           true,
                           filter_out_bf);
  mins_found = 0;
  while ((record = indexlr6.read())) {
    for (const auto& min : record.minimizers) {
      for (const auto h : filter_out_hashes) {
        TEST_ASSERT_NE(min.min_hash, h);
      }
      mins_found++;
    }
  }
  TEST_ASSERT_GE(mins_found, filter_in_hashes.size());

  btllib::Indexlr indexlr7(btllib::get_dirname(__FILE__) + "/indexlr.fq",
                           100,
                           5,
                           btllib::Indexlr::Flag::FILTER_IN |
                             btllib::Indexlr::Flag::FILTER_OUT |
                             btllib::Indexlr::Flag::SHORT_MODE,
                           3,
                           true,
                           filter_in_bf,
                           filter_out_bf);
  mins_found = 0;
  while ((record = indexlr7.read())) {
    for (const auto& min : record.minimizers) {
      bool found = false;
      for (const auto h : filter_in_hashes) {
        if (min.min_hash == h) {
          found = true;
          break;
        }
      }
      TEST_ASSERT(found);
      for (const auto h : filter_out_hashes) {
        TEST_ASSERT_NE(min.min_hash, h);
      }
      mins_found++;
    }
  }
  TEST_ASSERT_GE(mins_found, filter_in_hashes.size());

  return 0;
}
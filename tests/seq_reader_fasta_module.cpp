#include "btllib/seq_reader.hpp"
#include "btllib/util.hpp"
#include "helpers.hpp"

#include <iostream>
#include <memory>

int
main()
{
  const char* ids[] = { "asdf", "ghjk", "asdf", "ghjk", "wert", "yuio" };
  const char* seqs[] = { "ACTG",
                         "TGCA",
                         "ACDEFGHIKLMNOPQRSTUVWY",
                         "YWVUTSRQPONMLKIHGFEDCA",
                         "acdefghiklmnopqrstuvwy",
                         "ywvutsrqponmlkihgfedca" };

  for (int iteration = 0; iteration < 3; iteration++) {
    std::cerr << "Iteration " << iteration + 1 << std::endl;

    std::cerr << "Test small FASTA" << std::endl;

    btllib::SeqReader reader(btllib::get_dirname(__FILE__) +
                               "/input.fa.gz.bz2.xz",
                             btllib::SeqReader::Flag::SHORT_MODE);
    TEST_ASSERT_EQ(reader.get_format(), btllib::SeqReader::Format::FASTA)

    size_t i = 0;
    for (const auto record : reader) {
      TEST_ASSERT_EQ(record.id, ids[i]);
      TEST_ASSERT_EQ(record.seq, seqs[i]);
      TEST_ASSERT_EQ(record.qual, "");
      i++;
    }
    TEST_ASSERT_EQ(i, 2);
  }

  for (int iteration = 0; iteration < 3; iteration++) {
    std::cerr << "Iteration " << iteration + 1 << std::endl;

    std::cerr << "Test small FAA" << std::endl;

    btllib::SeqReader reader(btllib::get_dirname(__FILE__) + "/input.faa",
                             btllib::SeqReader::Flag::SHORT_MODE);
    TEST_ASSERT_EQ(reader.get_format(), btllib::SeqReader::Format::FASTA)

    size_t i = 2;
    for (const auto record : reader) {
      TEST_ASSERT_EQ(record.id, ids[i]);
      TEST_ASSERT_EQ(record.seq, seqs[i]);
      TEST_ASSERT_EQ(record.qual, "");
      i++;
    }
    TEST_ASSERT_EQ(i, 6);
  }

  return 0;
}
#include "btllib/randseq.hpp"

#include "helpers.hpp"

#include <iostream>
#include <regex>

int
main()
{
  {
    PRINT_TEST_NAME("DNA sequence generation")
    btllib::RandSeq rnd(btllib::RandSeq::SeqType::DNA,
                        btllib::RandSeq::Masking::NONE);
    TEST_ASSERT(std::regex_search(rnd.generate(100), std::regex("^[ACGT]*$")))
  }
  {
    PRINT_TEST_NAME("RNA sequence generation")
    btllib::RandSeq rnd(btllib::RandSeq::SeqType::RNA,
                        btllib::RandSeq::Masking::NONE);
    TEST_ASSERT(std::regex_search(rnd.generate(100), std::regex("^[ACGU]*$")))
  }
  {
    PRINT_TEST_NAME("random seed functionality")
    btllib::RandSeq rnd1(btllib::RandSeq::SeqType::DNA,
                         btllib::RandSeq::Masking::NONE);
    btllib::RandSeq rnd2(btllib::RandSeq::SeqType::DNA,
                         btllib::RandSeq::Masking::NONE);
    TEST_ASSERT(rnd1.generate(100) != rnd2.generate(100))
    rnd1.set_seed(42);
    rnd2.set_seed(42);
    TEST_ASSERT_EQ(rnd1.generate(100), rnd2.generate(100))
    rnd1.set_seed(1024);
    TEST_ASSERT(rnd1.generate(100) != rnd2.generate(100))
  }
}
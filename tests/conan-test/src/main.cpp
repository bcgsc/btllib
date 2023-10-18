#include <btllib/counting_bloom_filter.hpp>

int
main(int argc, char* argv[])
{
  btllib::KmerCountingBloomFilter8 cbf(1024, 3, 5);
  cbf.insert("ACGTA");
  cbf.insert("ACGTA");
  cbf.insert("ACGTA");
  cbf.insert("CCCCC");
  cbf.insert("TTTTT");
  return cbf.contains("ACGTA") == 3 ? EXIT_SUCCESS : EXIT_FAILURE;
}
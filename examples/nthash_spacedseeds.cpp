#include "../include/btllib/nthash.hpp"
#include "../include/btllib/seq_reader.hpp"

#include <iostream>

int
main()
{
  std::string seq = "ACGTACGT";
  std::vector<std::string> seeds = { "1001", "1111" };
  unsigned num_hashes = 2;
  btllib::SeedNtHash nthash(seq, seeds, num_hashes, 4);
  while (nthash.roll()) {
    for (int i = 0; i < seeds.size() * num_hashes; i++) {
      std::cout << std::hex << "0x" << nthash.hashes()[i] << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}
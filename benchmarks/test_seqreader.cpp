#include "btllib/seq_reader.hpp"

#include <cmath>
#include <iostream>
#include <omp.h>
#include <chrono>
#include <thread>

std::chrono::milliseconds WORK_TIME(2);

void work()
{
  std::this_thread::sleep_for(WORK_TIME);
}

int
main(int argc, char** argv)
{
  if (argc != 3) {
    std::cerr << "Missing args." << std::endl;
    std::exit(-1);
  }

  btllib::SeqReader reader(argv[1], btllib::SeqReader::Flag::LONG_MODE);
  omp_set_num_threads(std::atoi(argv[2]));

  int n = 0, slen = 0, qlen = 0;
#pragma omp parallel reduction(+ : n, slen, qlen)
  for (const auto record : reader) {
    work();
    ++n, slen += record.seq.size(), qlen += record.qual.size();
  }
  std::cerr << n << '\t' << slen << '\t' << qlen << std::endl;

  return 0;
}

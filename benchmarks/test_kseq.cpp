#include "kseq.h"
#include <stdio.h>
#include <zlib.h>
KSEQ_INIT(gzFile, gzread)

#include <cmath>
#include <iostream>
#include <omp.h>
#include <string>
#include <chrono>
#include <thread>

std::chrono::milliseconds WORK_TIME(5);

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

  gzFile fp;
  kseq_t* seq;
  fp = gzopen(argv[1], "r");
  seq = kseq_init(fp);

  omp_set_num_threads(std::atoi(argv[2]));

  int n = 0, slen = 0, qlen = 0;
#pragma omp parallel shared(fp, seq) reduction(+ : n, slen, qlen)
  {
    std::string seq_str;
    std::string qual_str;
    bool success = false;
    while (true) {
#pragma omp critical
      {
        success = (kseq_read(seq) >= 0);
        if (success) {
          seq_str = seq->seq.s;
          // qual_str = seq->qual.s;
        }
      }
      if (!success) {
        break;
      }
      work();
      ++n, slen += seq_str.size(), qlen += qual_str.size();
    }
  }
  std::cerr << n << '\t' << slen << '\t' << qlen << std::endl;

  kseq_destroy(seq);
  gzclose(fp);
  return 0;
}

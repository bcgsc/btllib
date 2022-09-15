#include "btllib/seq_reader.hpp"
#include "btllib/counting_bloom_filter.hpp"

#include <omp.h>

#include <iostream>
#include <vector>
#include <type_traits>
#include <string>

const size_t                  CBF_SIZE = 3ULL * 1024ULL * 1024ULL * 1024ULL;
const size_t                  SEEN_KMERS_BF_SIZE = 50ULL * 1024ULL * 1024ULL * 1024ULL;
const unsigned long           MAX_READS = 50000;
const unsigned                HASH_NUM = 4;
const unsigned                K = 25;
const std::vector<unsigned>   THREADS_NUMS = { 1, 16, 32, 64, 128 };
const size_t                  COUNTER_BYTES = 2;

template<typename T>
void fill_cbf(T& cbf, unsigned threads, const std::string& seqs_path) {
  std::cerr << "Filling..." << std::endl;

  omp_set_num_threads(threads);

  btllib::SeqReader reader(seqs_path, btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel shared(reader, cbf)
  for (auto record : reader) {
    if (record.num >= MAX_READS) { break; }
    if (record.seq.size() < K) { continue; }
    cbf.insert(record.seq);
  }
  std::cerr << "CBF FPRs = " << cbf.get_fpr(1) << ", " << cbf.get_fpr(2) << ", " << cbf.get_fpr(3) << std::endl;

  omp_set_num_threads(1);

  std::cerr << "Filled" << std::endl;
}

template<typename T>
void compare_cbfs(const T& cbf_singlethreaded, const T& cbf_multithreaded, unsigned threads, const std::string& seqs_path) {
  std::cerr << "Comparing..." << std::endl;

  btllib::check_error(cbf_singlethreaded.get_bytes() != cbf_multithreaded.get_bytes(), "Diff num of bytes!");

  omp_set_num_threads(threads);

  unsigned long total_counters = cbf_singlethreaded.get_bytes() / COUNTER_BYTES;
  unsigned long diff_counters = 0;
  unsigned long diff_sum = 0;

  btllib::KmerBloomFilter seen_kmers_bf(SEEN_KMERS_BF_SIZE, HASH_NUM, K);

  btllib::SeqReader reader(seqs_path, btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel shared(reader, cbf_singlethreaded, cbf_multithreaded) reduction(+:diff_counters, diff_sum)
  for (auto record : reader) {
    if (record.num >= MAX_READS) { break; }
    if (record.seq.size() < K) { continue; }
    btllib::NtHash nthash(record.seq, HASH_NUM, K);
    while (nthash.roll()) {
      bool seen = seen_kmers_bf.contains(nthash.hashes());
      if (seen) { continue; }

#pragma omp critical (seen_kmers_bf)
{
      seen = seen_kmers_bf.contains(nthash.hashes());
      if (!seen) {
        seen_kmers_bf.insert(nthash.hashes());
      }
}
      if (seen) { continue; }

      const uint64_t count_singlethreaded = cbf_singlethreaded.contains(nthash.hashes());
      const uint64_t count_multithreaded = cbf_multithreaded.contains(nthash.hashes());

      if (count_singlethreaded != count_multithreaded) {
        diff_counters++;

        const uint64_t diff = std::abs(long(count_singlethreaded) - long(count_multithreaded));
        diff_sum += diff;
      }
    }
  }
  std::cerr << "Seen kmers BF FPR = " << seen_kmers_bf.get_fpr() << std::endl;

  std::cerr << "Total counters = " << total_counters << std::endl;
  std::cerr << "Diff counters = " << diff_counters << std::endl;
  std::cerr << "Diff sum = " << diff_sum << std::endl;
  std::cerr << "Avg counter diff = " << (double(diff_sum) / double(total_counters)) << std::endl;

  omp_set_num_threads(1);

  std::cerr << "Compared" << std::endl;
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Usage: test_threadsafety <threads> <seqs>" << std::endl;
    return EXIT_FAILURE;
  }
  const unsigned threads = std::stoul(argv[1]);
  const std::string seqs_path = argv[2];

  btllib::KmerCountingBloomFilter16 cbf_singlethreaded(CBF_SIZE, HASH_NUM, K);
  fill_cbf(cbf_singlethreaded, 1, seqs_path);

  std::cerr << "Running multithreaded scenarios...\n"
            << "Threads = " << threads << std::endl;

  btllib::KmerCountingBloomFilter16 cbf_multithreaded(CBF_SIZE, HASH_NUM, K);
  fill_cbf(cbf_multithreaded, threads, seqs_path);
  compare_cbfs(cbf_singlethreaded, cbf_multithreaded, 30, seqs_path);

  return 0;
}

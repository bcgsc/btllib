#include "btllib/mi_bloom_filter.hpp"
#include "btllib/nthash.hpp" //parse_seeds
#include "btllib/seq_reader.hpp"
#include "btllib/status.hpp"
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <omp.h>

#include <limits>
#include <string>
#include <vector>

const static std::string PROGNAME = "mibf";
const static std::string VERSION = "1.0";
const static size_t DEFAULT_THREADS = 5;
const static size_t DEFAULT_SIZE = 1000000000;
const static double DEFAULT_OCCUPANCY = 0.5;
const static size_t DEFAULT_SEQ_READER_THREADS = 6;

using ID_type = uint16_t;
using SpacedSeed = std::vector<unsigned>;

static void
print_usage()
{
  std::cerr
    << "Usage: " << PROGNAME
    << " -p MI_BF_PREFIX [OPTION]... [FILE]\n"
       "  -p file_prefix     filter prefix and filter ID, required option.\n"
       "--------------------------------------------------------------\n"
       "  -k kmer_size       k-mer size.\n"
       "  -g hash_num        hash number.\n"
       "  or\n"
       "  -s spaced_seeds    expects list of seed 1s & 0s separated by "
       "spaces.\n"
       "                     required if '-k' and '-g' are not given.\n"
       "--------------------------------------------------------------\n"
       "  -m mi_bf_size      output mi_bf bit vector size in bits, default is "
       "10^9.\n"
       "  or\n"
       "  -b occupancy       occupancy of Bloom filter, default is 0.5.\n"
       "                     required if '-m' is not given.\n"
       "  -n number_of_elems the number of expected distinct k-mer frames,\n"
       "by default determined from sequences within files.\n"
       "                     required if '-m' is not given.\n"
       "--------------------------------------------------------------\n"
       "  -f by_file         assign IDs by file rather than by fasta header.\n"
       "  -t threads         number of threads (default 5, max 32)\n"
       "  -v verbose         show verbose output.\n"
       "  --help             display this help and exit.\n"
       "  --version          display version and exit.\n"
       "  FILE               space separated list of FASTA/Q files."
    << std::endl;
}

static std::vector<std::string>
split_string_by_space(const std::string& input_string)
{
  std::vector<std::string> current_string;
  std::string temp;
  std::stringstream converter(input_string);
  while (converter >> temp) {
    current_string.push_back(temp);
  }
  assert(!current_string.empty());
  return current_string;
}

static unsigned
assert_id_size_and_count_kmers(const std::vector<std::string>& read_paths,
                               const bool& by_file,
                               const unsigned& kmer_size)
{
  unsigned total_id = 0;
  unsigned expected_elements = 0;

  for (const auto& read_path : read_paths) {
    btllib::SeqReader reader(read_path, btllib::SeqReader::Flag::SHORT_MODE);
    for (const auto& record : reader) {
      if (!by_file && !record.seq.empty()) {
        total_id++;
      }
    }
  }

  for (const auto& read_path : read_paths) {
    btllib::SeqReader reader(read_path, btllib::SeqReader::Flag::SHORT_MODE);
    if (by_file) {
      total_id++;
    }
    for (const auto& record : reader) {
      if (!by_file) {
        total_id++;
      }
      expected_elements +=
        record.seq.size() > kmer_size ? record.seq.size() - kmer_size : 0;
    }
  }
  btllib::check_error(
    !(total_id < static_cast<unsigned>(std::numeric_limits<ID_type>::max())),
    "Total ID number overflows ID type.");
  return expected_elements;
}

template<typename T>
inline void
insert_to_bv(btllib::MIBloomFilter<ID_type>& mi_bf,
             T& nthash,
             int& mi_bf_stage,
             ID_type& id)
{
  while (nthash.roll()) {
    if (mi_bf_stage == 0) {
      mi_bf.insert_bv(nthash.hashes());
    } else if (mi_bf_stage == 1) {
      mi_bf.insert_id(nthash.hashes(), id);
    } else if (mi_bf_stage == 2) {
      mi_bf.insert_saturation(nthash.hashes(), id);
    }
  }
}

int
main(int argc, char* argv[])
{
  try {
    int c;
    bool failed = false;
    int optindex = 0;
    // static int version = 0;
    static int help = 0, version = 0;
    bool verbose = false;
    int thread_count = DEFAULT_THREADS;
    double occupancy = DEFAULT_OCCUPANCY;
    unsigned mi_bf_size = DEFAULT_SIZE, kmer_size = 0, hash_num = 0,
             expected_elements = 0;

    bool spaced_seed_set = false;
    bool output_path_set = false;
    bool by_file = false;
    bool occupancy_set = false;

    std::string output_path;

    std::vector<std::string> spaced_seeds_string;
    std::vector<SpacedSeed> spaced_seeds;

    static const struct option longopts[] = {
      { "help", no_argument, &help, 1 },
      { "version", no_argument, &version, 1 },
      { nullptr, 0, nullptr, 0 }
    };

    while ((c = getopt_long(argc, // NOLINT(concurrency-mt-unsafe)
                            argv,
                            "p:k:g:s:m:n:b:ft:v",
                            longopts,
                            &optindex)) != -1) {
      switch (c) {
        case 0:
          break;
        case 'p':
          output_path = optarg;
          output_path_set = true;
          break;
        case 'k':
          kmer_size = std::stoul(optarg);
          break;
        case 'g':
          hash_num = std::stoul(optarg);
          break;
        case 's':
          spaced_seed_set = true;
          spaced_seeds_string = split_string_by_space(optarg);
          spaced_seeds = btllib::parse_seeds(spaced_seeds_string);
          hash_num = spaced_seeds.size();
          kmer_size = spaced_seeds_string[0].size();
          break;
        case 'b':
          occupancy = std::stod(optarg);
          occupancy_set = true;
          break;
        case 'f':
          by_file = true;
          break;
        case 'm':
          mi_bf_size = std::stoul(optarg);
          break;
        case 'n':
          expected_elements = std::stoul(optarg);
          break;
        case 't':
          thread_count = std::stoi(optarg);
          break;
        case 'v':
          verbose = true;
          break;
        default:
          std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
      }
    }
    if (kmer_size == 0) {
      btllib::log_error("missing option -- 'k'");
      failed = true;
    }
    if (!failed && hash_num == 0) {
      btllib::log_error("missing option -- 'h'");
      failed = true;
    }
    if (!failed && !output_path_set) {
      btllib::log_error("missing option -- 'p'");
      failed = true;
    }
    if (verbose) {
      std::cout << "verbose" << std::endl;
    }

    std::vector<std::string> read_paths(&argv[optind], &argv[argc]);
    if (argc < 2 || failed) {
      std::cout << std::endl;
      print_usage();
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
    std::map<std::string, ID_type> read_ids;
    unsigned kmer_count =
      assert_id_size_and_count_kmers(read_paths, by_file, kmer_size);

    if (expected_elements > 0 || occupancy_set) {
      mi_bf_size = btllib::MIBloomFilter<ID_type>::calc_optimal_size(
        kmer_count, hash_num, occupancy);
      btllib::log_info("Optimal size is calculated: " +
                       std::to_string(mi_bf_size));
    }

    // set thread number configuration
#if defined(_OPENMP)
    if (thread_count > 0) {
      omp_set_num_threads(thread_count);
    }
#endif

    btllib::MIBloomFilter<ID_type> mi_bf(mi_bf_size, hash_num);
    const char* stages[3] = { "BV Insertion", "ID Insertion", "Saturation" };
    ID_type id_counter = 0;
    std::map<std::string, ID_type> ids;
    for (int mi_bf_stage = 0; mi_bf_stage < 3; mi_bf_stage++) {
      btllib::log_info(stages[mi_bf_stage] + std::string(" stage started"));

      for (auto& read_path : read_paths) {
        btllib::SeqReader reader(read_path,
                                 btllib::SeqReader::Flag::SHORT_MODE,
                                 DEFAULT_SEQ_READER_THREADS);
#pragma omp parallel default(none) shared(ids,                                 \
                                          id_counter,                          \
                                          mi_bf,                               \
                                          mi_bf_stage,                         \
                                          hash_num,                            \
                                          kmer_size,                           \
                                          by_file,                             \
                                          spaced_seed_set,                     \
                                          spaced_seeds,                        \
                                          reader)
        try {
          for (const auto record : reader) {
#pragma omp critical
            {
              if (ids.find(record.id) == ids.end()) {
                ids[record.id] = !by_file ? id_counter++ : id_counter;
              }
            }
            if (spaced_seed_set) {
              btllib::SeedNtHash nthash(record.seq, spaced_seeds, 1, kmer_size);
              insert_to_bv<btllib::SeedNtHash>(
                mi_bf, nthash, mi_bf_stage, ids[record.id]);
            } else {
              btllib::NtHash nthash(record.seq, hash_num, kmer_size);
              insert_to_bv<btllib::NtHash>(
                mi_bf, nthash, mi_bf_stage, ids[record.id]);
            }
          }
          if (by_file) {
            id_counter++;
          }
        } catch (const std::exception& e) {
          btllib::log_error(e.what());
          std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
        }
      }
      if (mi_bf_stage == 0) {
        mi_bf.complete_bv_insertion();
      }
      if (mi_bf_stage == 1) {
        mi_bf.complete_id_insertion();
      }
    }
    mi_bf.save(output_path);

    std::ofstream ids_file;
    ids_file.open(output_path + ".ids");
    for (auto& elem : ids) {
      ids_file << elem.first << "\t" << elem.second << std::endl;
    }
    ids_file.close();
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
}

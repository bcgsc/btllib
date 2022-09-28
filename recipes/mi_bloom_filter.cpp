#include "btllib/mi_bloom_filter.hpp"
#include "btllib/status.hpp"
#include "btllib/seq_reader.hpp"

#include <getopt.h>

#include <string>
#include <vector>

const static std::string PROGNAME = "mi_bloom_filter";
const static std::string VERSION = "1.0";
const static size_t MAX_THREADS = 5;
const static size_t DEFAULT_THREADS = MAX_THREADS;

static void
print_error_msg(const std::string& msg)
{
  std::cerr << PROGNAME << ' ' << VERSION << ": " << msg << std::endl;
}

static void
print_usage()
{
  std::cerr << "Usage: " << PROGNAME
	<< "-k K -h H -o output"
	<< "[-i id_file] [-m mi_bf_size] FILE...\n\n"
		"  -k K        k-mer size.\n"
		"  -h H        hash number.\n"
		"  -o output   Path of mi_bf."
		"  -i id_file  Get IDs from file. Default is by file order."
		"  -m mi_bf_size  Output mi_bf bit vector size in bits. Default is 10^9."
		"  -t T        Use T number of threads (default 5, max 5) per "
		"  -v          Show verbose output.\n"
		"  --help      Display this help and exit.\n"
		"  --version   Display version and exit.\n"
		"  FILE        Space separated list of FASTA/Q files."    
	<< std::endl;
}

template<typename T>
inline void insert_to_bv(btllib::MIBloomFilter<uint8_t>& mi_bf, T& nthash, int& mi_bf_stage, uint8_t& ID){
	while (nthash.roll()) {
		if(mi_bf_stage == 0){
			mi_bf.insert_bv(nthash.hashes());
		} else if(mi_bf_stage == 1){
			mi_bf.insert_id(nthash.hashes(), ID);
		} else if(mi_bf_stage == 2){
			mi_bf.insert_saturation(nthash.hashes(), ID);
        	}
	}
}

int
main(int argc, char* argv[])
{
  try {
    int c;
    int optindex = 0;
    int help = 0, version = 0;
    bool verbose = false;
    int mi_bf_size = 1000000000;

    unsigned k = 0, h = 0, t = DEFAULT_THREADS;
    bool k_set = false;
    bool h_set = false;
    bool output_bf_set = false;
    bool spaced_seed_set = false;
    std::string output_path, id_file;
    
    std::vector<std::string> spaced_seeds_string;
    std::vector<btllib::SpacedSeed> spaced_seeds;

    static const struct option longopts[] = {
      { "help", no_argument, &help, 1 },
      { "version", no_argument, &version, 1 },
      { nullptr, 0, nullptr, 0 }
    };

      while ((c = getopt_long(argc, // NOLINT(concurrency-mt-unsafe)
                            argv,
                            "k:h:o:i:m:t:v",
                            longopts,
                            &optindex)) != -1) {
      switch (c) {
        case 0:
          break;
        case 'k':
          k_set = true;
          k = std::stoul(optarg);
          break;
        case 'h':
          h_set = true;
          h = std::stoul(optarg);
          break;
        case 'o':
          output_path = optarg;
          output_bf_set = true;
	  break;
        case 'i':
          id_file = optarg;
          break;
        case 'm':
          mi_bf_size = std::stoul(optarg);
          break;
        case 't':
          t = std::stoul(optarg);
          break;
        case 'v':
          verbose = true;
          break;
        default:
          std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
      }
    }
    if (!k_set) {
      print_error_msg("missing option -- 'k'");
    //  failed = true;
    }
    if (!output_bf_set) {
      print_error_msg("missing option -- 'o'");
    //  failed = true;
    }
    if (!h_set) {
      print_error_msg("missing option -- 'h'");
    //  failed = true;
    }
    if(verbose){std::cout << "verbose" << std::endl;}
    t = t ? t : t;
    std::vector<std::string> read_paths(&argv[optind], &argv[argc]);
    if (argc < 2) {
      print_usage();
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }

    btllib::MIBloomFilter<uint8_t> mi_bf(mi_bf_size, h);

    for(int mi_bf_stage = 0; mi_bf_stage < 2; mi_bf_stage++){
		if(mi_bf_stage == 0){std::cout << "BV Insertion stage started" << std::endl;}
		else if(mi_bf_stage == 1){std::cout << "ID Insertion stage" << std::endl;}
		else if(mi_bf_stage == 2){std::cout << "Saturation stage started" << std::endl;}

    		for(uint m = 1; m <= read_paths.size(); m++){
			uint8_t ID = m;
			
			btllib::SeqReader reader(read_paths[m], btllib::SeqReader::Flag::SHORT_MODE, 6);

			for (const auto record : reader) {
				if(spaced_seed_set){
					btllib::SeedNtHash nthash(record.seq, spaced_seeds, 1, spaced_seeds_string[0].size());
        	       	                insert_to_bv<btllib::SeedNtHash>(mi_bf, nthash, mi_bf_stage, ID);
				} else {
					btllib::NtHash nthash(record.seq, h, k);
        	       	                insert_to_bv<btllib::NtHash>(mi_bf, nthash, mi_bf_stage, ID);
				}
			}
		}
		if(mi_bf_stage == 0){ mi_bf.complete_bv_insertion();}
   }
	mi_bf.save(output_path);
   }    
   catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
   }

}

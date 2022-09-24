#include "btllib/mi_bloom_filter.hpp"
#include "btllib/status.hpp"

#include <string>
#include <vector>

const static std::string PROGNAME = "mi_bloom_filter";
const static std::string VERSION = btllib::PROJECT_VERSION;

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
	<< "[-i id_file] FILE...\n\n"
		"  -k K        k-mer size.\n"
		"  -h H        hash number.\n"
		"  -o output   Path of mi_bf."
		"  -i id_file  Get IDs from file. Default is by file order."
		"  -t T        Use T number of threads (default 5, max 5) per "
		"  -v          Show verbose output.\n"
		"  --help      Display this help and exit.\n"
		"  --version   Display version and exit.\n"
		"  FILE        Space separated list of FASTA/Q files."    
	<< std::endl;
}

template<typename T>
inline void insert_to_bv(btllib::MIBloomFilter<uint8_t>& mi_bf, T& nthash, int& mi_bf_stage, uint8_t& ID, size_t& hash_num){
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
    int long_mode = 0;

    unsigned k = 0, h = 0, t = DEFAULT_THREADS;
    bool k_set = false;
    bool h_set = false;
    bool output_bf_set = false;
    std::string output_path;
    static const struct option longopts[] = {
      { "help", no_argument, &help, 1 },
      { "version", no_argument, &version, 1 },
      { nullptr, 0, nullptr, 0 }
    };

while ((c = getopt_long(argc, // NOLINT(concurrency-mt-unsafe)
                            argv,
                            "k:h:o:i:t:v",
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
    std::vector<std::string> infiles(&argv[optind], &argv[argc]);
    if (argc < 2) {
      print_usage();
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }

    btllib::MIBloomFilter<uint8_t> mi_bf(mi_bf_size, h);

    for(int mi_bf_stage = 0; mi_bf_stage < stage_count; mi_bf_stage++){
		if(mi_bf_stage == 0){std::cout << "BV Insertion stage started" << std::endl;}
		else if(mi_bf_stage == 1){std::cout << "ID Insertion stage" << std::endl;}
		else if(mi_bf_stage == 2){std::cout << "Saturation stage started" << std::endl;}
    
    		for(int m = 1; m <= read_paths.size(); m++){
			btllib::SeqReader reader(read, btllib::SeqReader::Flag::SHORT_MODE, 6);
			for (const auto record : reader) {
				if(params.spaced_seed_bool){
					btllib::SeedNtHash nthash(record.seq, spaced_seeds, 1, spaced_seeds_string[0].size());
        	       	                insert_to_bv<btllib::SeedNtHash>(mi_bf, nthash, mi_bf_stage, ID, params.hash_number);
				} else {
					btllib::NtHash nthash(record.seq, params.hash_number, params.kmer_size);
        	       	                insert_to_bv<btllib::NtHash>(mi_bf, nthash, mi_bf_stage, ID, params.hash_number);
				}
			}
		}
		if(mi_bf_stage == 0){ mi_bf.complete_bv_insertion();}
	}
	mi_bf.save(output);
  }
}










#include "btllib/randseq.hpp"
#include "btllib/seq_writer.hpp"
#include "btllib/status.hpp"
#include "config.hpp"

#include <argparse/argparse.hpp>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <random>
#include <thread>

struct Arguments
{
  btllib::RandSeq::SeqType seq_type;
  btllib::RandSeq::Masking mask;
  unsigned num_sequences;
  int num_threads;
  size_t min_length, max_length;
  std::string out_path;

  Arguments(int argc, char** argv)
  {
    argparse::ArgumentParser parser("randseq", btllib::PROJECT_VERSION);

    parser.add_argument("-s")
      .help("Sequence type (dna | rna)")
      .default_value(std::string("dna"));

    parser.add_argument("-m")
      .help("Masking type (none | soft | hard)")
      .default_value(std::string("none"));

    parser.add_argument("-n")
      .help("Number of sequences to generate")
      .default_value(1U)
      .scan<'u', unsigned>();

    parser.add_argument("-l")
      .help(
        "Sequence length. To generate sequences with random lengths from the "
        "range [a, b], use a:b")
      .required();

    parser.add_argument("-t")
      .help("Number of parallel threads")
      .default_value(1)
      .scan<'i', int>();

    parser.add_argument("-o")
      .help("Path to output file. Use '-' to write to stdout.")
      .default_value(std::string("-"));

    try {
      parser.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
      std::cerr << err.what() << std::endl;
      std::cerr << parser;
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }

    auto seq_type_str = parser.get("-s");
    if (seq_type_str == "dna") {
      seq_type = btllib::RandSeq::SeqType::DNA;
    } else if (seq_type_str == "rna") {
      seq_type = btllib::RandSeq::SeqType::RNA;
    } else {
      btllib::check_error(true,
                          "Invalid sequence type: " + seq_type_str +
                            ". Should be 'dna' or 'rna'.");
    }

    auto mask_str = parser.get("-m");
    if (mask_str == "none") {
      mask = btllib::RandSeq::Masking::NONE;
    } else if (mask_str == "soft") {
      mask = btllib::RandSeq::Masking::SOFT;
    } else if (mask_str == "hard") {
      mask = btllib::RandSeq::Masking::HARD;
    } else {
      btllib::check_error(true,
                          "Invalid masking option: " + mask_str +
                            ". Should be 'none', 'soft', or 'hard'.");
    }

    auto length_range = parser.get("-l");
    auto i_sep = length_range.find(':');
    if (i_sep == std::string::npos) {
      min_length = std::stoull(length_range);
      max_length = std::stoull(length_range);
    } else {
      auto a = length_range.substr(0, i_sep);
      auto b = length_range.substr(i_sep + 1, length_range.size() - i_sep - 1);
      min_length = std::stoull(a);
      max_length = std::stoull(b);
    }

    num_sequences = parser.get<unsigned>("-n");
    num_threads = parser.get<int>("-t");
    out_path = parser.get("-o");
  }
};

class Logger
{
  static constexpr double LOG_DELAY_SECONDS = 0.1;
  std::chrono::time_point<std::chrono::system_clock> start_time, last_log;
  uint64_t generated_bp = 0, last_generated_bp = 0;
  unsigned generated_seqs = 0, total_seqs;
  size_t last_line_length = 0;
  const std::string& unit;

public:
  Logger(unsigned num_sequences, const std::string& unit)
    : start_time(std::chrono::system_clock::now())
    , last_log(std::chrono::system_clock::now())
    , total_seqs(num_sequences)
    , unit(unit)
  {
  }

  void add_sequence(size_t seq_len)
  {
    generated_bp += seq_len;
    ++generated_seqs;
    auto current_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = (current_time - last_log);
    if (elapsed.count() > LOG_DELAY_SECONDS) {
      long double diff_bp = generated_bp - last_generated_bp;
      auto speed = std::llround(diff_bp / elapsed.count());
      auto progress =
        std::lround((double)generated_seqs / (double)total_seqs * 100.0);
      double seq_ratio =
        ((double)total_seqs - (double)generated_seqs) / (double)generated_seqs;
      std::chrono::duration<double> total_elapsed = (current_time - start_time);
      auto remaining = std::lround(seq_ratio * total_elapsed.count() + 1);
      std::string line = "[" + std::to_string(generated_seqs) + "/" +
                         std::to_string(total_seqs) + "] Generated " +
                         std::to_string(generated_bp) + unit + " @" +
                         std::to_string(speed) + unit + "/s (" +
                         std::to_string(progress) + "%, ~" +
                         std::to_string(remaining) + "s remaining)";
      std::cerr << line << "\r";
      last_log = current_time;
      last_line_length = line.size();
      last_generated_bp = generated_bp;
    }
  }

  void stop()
  {
    auto end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = (end_time - start_time);
    std::cerr << std::string(last_line_length, ' ') << "\r";
    std::cerr << "Generated " << generated_seqs << " sequence(s) (total "
              << generated_bp << unit << ") in " << elapsed.count() << "s"
              << std::endl;
  }
};

int
main(int argc, char** argv)
{
  try {
    Arguments args(argc, argv);
    omp_set_num_threads(args.num_threads);

    std::random_device rd;
    std::default_random_engine rng(rd());
    std::uniform_int_distribution<size_t> dist(args.min_length,
                                               args.max_length);
    btllib::SeqWriter writer(args.out_path);
    btllib::RandSeq rnd(args.seq_type, args.mask);

    std::string unit =
      args.seq_type == btllib::RandSeq::SeqType::PROTEIN
        ? "aa"
        : "bp";
    Logger log(args.num_sequences, unit);

#pragma omp parallel for shared(args, rnd, rng, dist, writer, log) default(none)
    for (unsigned i = 0; i < args.num_sequences; i++) {
      std::string seq = rnd.generate(dist(rng));
      writer.write(std::to_string(i + 1), "", seq, "");
#pragma omp critical
      log.add_sequence(seq.size());
    }

    log.stop();
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }

  return 0;
}
#include "../include/btllib/rolling_hash.hpp"

#include <cassert>
#include <iostream>
#include <string>

int
main()
{
  std::string kmer = "ACGTACACTGGACTGAGTCT";

  {
    std::cerr << "Invariant hash values" << std::endl;
    btllib::RollingHash rolling_hash(kmer, kmer.size(), 3);

    /* Hash values*/
    const std::vector<uint64_t> hashes = { 10434435546371013747U,
                                           16073887395445158014U,
                                           8061578976118370557U };

    rolling_hash.roll();
    assert(rolling_hash.get_hash_num() == hashes.size());
    size_t i;
    for (i = 0; i < rolling_hash.get_hash_num(); ++i) {
      assert(rolling_hash.hashes()[i] == hashes[i]);
    }
    assert(i == 3);
  }

  std::cerr << "Reverse complement" << std::endl;
  {
    /* Reverse complement of kmer*/
    std::string rc_kmer = "AGACTCAGTCCAGTGTACGT";

    btllib::RollingHash rolling_hash(kmer, 20, 3);
    btllib::RollingHash rolling_hash_rc(rc_kmer, 20, 3);

    rolling_hash.roll();
    rolling_hash_rc.roll();
    assert(rolling_hash.get_hash_num() == rolling_hash_rc.get_hash_num());
    size_t i;
    for (i = 0; i < rolling_hash.get_hash_num(); ++i) {
      assert(rolling_hash.hashes()[i] == rolling_hash_rc.hashes()[i]);
    }
    assert(i == 3);
  }

  std::cerr << "Rolling hash values" << std::endl;
  {
    btllib::RollingHash rolling_hash(kmer, 18, 3);

    /* 18-mers of kmer*/
    std::string kmer1 = "ACGTACACTGGACTGAGT";
    std::string kmer2 = "CGTACACTGGACTGAGTC";
    std::string kmer3 = "GTACACTGGACTGAGTCT";

    std::vector<btllib::RollingHash> rolling_hash_vector = {
      btllib::RollingHash(kmer1, kmer1.size(), rolling_hash.get_hash_num()),
      btllib::RollingHash(kmer2, kmer2.size(), rolling_hash.get_hash_num()),
      btllib::RollingHash(kmer3, kmer3.size(), rolling_hash.get_hash_num())
    };

    size_t i;
    for (i = 0; rolling_hash.roll() && rolling_hash_vector[i].roll(); ++i) {
      for (size_t j = 0; j < rolling_hash.get_hash_num(); ++j) {
        assert(rolling_hash.hashes()[j] == rolling_hash_vector[i].hashes()[j]);
      }
    }
    assert(i == 3);
  }

  std::cerr << "Spaced seeds" << std::endl;
  {
    std::vector<std::string> seeds = { "11111100000000111111",
                                       "11111111000011111111" };

    /* Point Mutations of Kmer*/
    std::string kmerM1 = "ACGTACACTTGACTGAGTCT";
    std::string kmerM2 = "ACGTACACTGTACTGAGTCT";
    std::string kmerM3 = "ACGTACACTGCACTGAGTCT";
    assert(kmerM1.size() == seeds[0].size());
    assert(kmerM1.size() == seeds[1].size());

    btllib::SeedRollingHash seed_rolling_hash(kmer, kmer.size(), seeds, 2);

    std::vector<btllib::SeedRollingHash> seed_rolling_hash_vector = {
      btllib::SeedRollingHash(kmerM1,
                              seed_rolling_hash.get_k(),
                              seeds,
                              seed_rolling_hash.get_hash_num_per_seed()),
      btllib::SeedRollingHash(kmerM2,
                              seed_rolling_hash.get_k(),
                              seeds,
                              seed_rolling_hash.get_hash_num_per_seed()),
      btllib::SeedRollingHash(kmerM3,
                              seed_rolling_hash.get_k(),
                              seeds,
                              seed_rolling_hash.get_hash_num_per_seed())
    };
    assert(seed_rolling_hash.get_hash_num() == seeds.size() * 2);
    assert(seed_rolling_hash.get_hash_num() ==
           seed_rolling_hash_vector[0].get_hash_num());

    seed_rolling_hash.roll();
    size_t i;
    for (i = 0; i < seed_rolling_hash_vector.size() &&
                seed_rolling_hash_vector[i].roll();
         i++) {
      for (size_t j = 0; j < seed_rolling_hash.get_hash_num(); j++) {
        assert(seed_rolling_hash.hashes()[j] ==
               seed_rolling_hash_vector[i].hashes()[j]);
      }
    }
    assert(i == 3);
  }

  std::cerr << "RNA" << std::endl;
  {
    btllib::RollingHash dna_rolling_hash(kmer, 20, 3);

    std::string rna_kmer = "ACGUACACUGGACUGAGUCU";
    btllib::RollingHash rna_rolling_hash(kmer, 20, 3);

    dna_rolling_hash.roll();
    rna_rolling_hash.roll();
    size_t i;
    for (i = 0; i < dna_rolling_hash.get_hash_num(); ++i) {
      assert(dna_rolling_hash.hashes()[i] == rna_rolling_hash.hashes()[i]);
    }
    assert(i == 3);
  }

  return 0;
}
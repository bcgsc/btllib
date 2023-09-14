#include "btllib/phred_nthash.hpp"
#include "btllib/seq.hpp"
#include "btllib/util.hpp"
#include "helpers.hpp"

#include <iostream>
#include <queue>
#include <stack>
#include <string>

int
main()
{

  {
    PRINT_TEST_NAME("Test skip due to base with low Phred score")

    std::string seq = "ACATGCCATGC";
    std::string qual = "11111011111";
    const unsigned k = 5;
    const unsigned h = 3;

    const std::vector<std::array<uint64_t, h>> hashes = {
      { 0xf59ecb45f0e22b9c, 0x4969c33ac240c129, 0x688d616f0d7e08c3 },
      { 0x38cc00f940aebdae, 0xab7e1b110e086fc6, 0x11a1818bcfdd553 }
    };

    btllib::PhredNtHash phred_nthash(seq, h, k, 16, qual);

    for (const auto& h_vals : hashes) {
      phred_nthash.roll();
      TEST_ASSERT_ARRAY_EQ(h_vals, phred_nthash.hashes(), h);
    }

    phred_nthash.roll_back();
    TEST_ASSERT_EQ(phred_nthash.get_pos(), 0);
    TEST_ASSERT_ARRAY_EQ(hashes[0], phred_nthash.hashes(), h);

  }

  {
    PRINT_TEST_NAME("k-mer hash values")

    std::string seq = "ACATGCATGCA";
    std::string qual = "$$%%)*0)'%%";
    const unsigned k = 5;
    const unsigned h = 3;

    const std::vector<std::array<uint64_t, h>> hashes = {
      { 0xf59ecb45f0e22b9c, 0x4969c33ac240c129, 0x688d616f0d7e08c3 },
      { 0x38cc00f940aebdae, 0xab7e1b110e086fc6, 0x11a1818bcfdd553 },
      { 0x603a48c5a11c794a, 0xe66016e61816b9c4, 0xc5b13cb146996ffe }
    };

    btllib::PhredNtHash phred_nthash(seq, h, k, 0, qual);

    for (const auto& h_vals : hashes) {
      phred_nthash.roll();
      TEST_ASSERT_ARRAY_EQ(h_vals, phred_nthash.hashes(), h);
    }
  }

  {
    PRINT_TEST_NAME("k-mer rolling")

    std::string seq = "AGTCAGTC";
    std::string qual = "$$%%)*0)";
    unsigned h = 3;
    unsigned k = 4;

    btllib::PhredNtHash phred_nthash(seq, h, k, 0, qual);
    std::vector<uint64_t*> hashes;

    while (phred_nthash.roll()) {
      uint64_t* h_vals = new uint64_t[h];
      std::copy(phred_nthash.hashes(), phred_nthash.hashes() + h, h_vals);
      hashes.push_back(h_vals);
    }

    TEST_ASSERT_EQ(hashes.size(), seq.length() - k + 1);

    // check same hash value for identical k-mers (first and last)
    TEST_ASSERT_ARRAY_EQ(hashes[0], hashes[hashes.size() - 1], h);
  }

  {
    PRINT_TEST_NAME("k-mer rolling vs ntbase hash values")

    std::string seq = "ACGTACACTGGACTGAGTCT";
    std::string qual = "$$%%)*0)'%%%%)*0)'%%";

    btllib::PhredNtHash phred_nthash(seq, 3, seq.size() - 2, 0, qual);
    /* 18-mers of kmer*/
    std::string kmer1 = "ACGTACACTGGACTGAGT";
    std::string kmer2 = "CGTACACTGGACTGAGTC";
    std::string kmer3 = "GTACACTGGACTGAGTCT";

    btllib::PhredNtHash phred_nthash_vector[] = {
      btllib::PhredNtHash(kmer1, phred_nthash.get_hash_num(), kmer1.size(), 0, qual),
      btllib::PhredNtHash(kmer2, phred_nthash.get_hash_num(), kmer2.size(), 0, qual),
      btllib::PhredNtHash(kmer3, phred_nthash.get_hash_num(), kmer3.size(), 0, qual)
    };

    size_t i;
    for (i = 0; phred_nthash.roll() && phred_nthash_vector[i].roll(); ++i) {
      for (size_t j = 0; j < phred_nthash.get_hash_num(); ++j) {
        TEST_ASSERT_EQ(phred_nthash.hashes()[j], phred_nthash_vector[i].hashes()[j]);
      }
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    PRINT_TEST_NAME("canonical hashing")

    std::string seq_f = "ACGTACACTGGACTGAGTCT";
    std::string seq_r = "AGACTCAGTCCAGTGTACGT";
    std::string qual = "$$%%)*0)'%%%%)*0)'%%";
    
    unsigned h = 3;

    btllib::PhredNtHash phred_nthash_f(seq_f, h, seq_f.size(), 0, qual);
    btllib::PhredNtHash phred_nthash_r(seq_r, h, seq_r.size(), 0, qual);

    phred_nthash_f.roll();
    phred_nthash_r.roll();
    TEST_ASSERT_EQ(phred_nthash_f.get_hash_num(), phred_nthash_r.get_hash_num())
    TEST_ASSERT_ARRAY_EQ(phred_nthash_f.hashes(), phred_nthash_r.hashes(), h)
  }

  {
    PRINT_TEST_NAME("k-mer back rolling")

    std::string seq = "ACTAGCTG";
    std::string qual = "$$%%)*0%";
    unsigned h = 3;
    unsigned k = 5;

    btllib::PhredNtHash phred_nthash(seq, h, k, 0, qual);
    std::stack<uint64_t*> hashes;

    while (phred_nthash.roll()) {
      uint64_t* h_vals = new uint64_t[h];
      std::copy(phred_nthash.hashes(), phred_nthash.hashes() + h, h_vals);
      hashes.push(h_vals);
    }

    TEST_ASSERT_EQ(hashes.size(), seq.length() - k + 1)

    do {
      TEST_ASSERT_ARRAY_EQ(phred_nthash.hashes(), hashes.top(), h)
      hashes.pop();
    } while (phred_nthash.roll_back());
  }

  {
    PRINT_TEST_NAME("skipping Ns")

    std::string seq = "ACGTACACTGGACTGAGTCT";
    std::string qual = "$$%%)*0)'%%%%)*0)'%%";
    std::string seq_with_ns = seq;

    TEST_ASSERT_GE(seq_with_ns.size(), 10)
    seq_with_ns[seq_with_ns.size() / 2] = 'N';
    seq_with_ns[seq_with_ns.size() / 2 + 1] = 'N';
    unsigned k = (seq.size() - 2) / 2 - 1;
    btllib::PhredNtHash phred_nthash(seq_with_ns, 3, k, 0, qual);

    std::vector<uint64_t> positions;
    for (size_t i = 0; i < seq_with_ns.size() / 2 - k + 1; i++) {
      positions.push_back(i);
    }
    for (size_t i = seq_with_ns.size() / 2 + 2; i < seq_with_ns.size() - k + 1;
         i++) {
      positions.push_back(i);
    }

    size_t i = 0;
    while (phred_nthash.roll()) {
      TEST_ASSERT_EQ(phred_nthash.get_pos(), positions[i])
      i++;
    }
    TEST_ASSERT_EQ(positions.size(), i)
  }

  {
    std::cerr << "Testing RNA" << std::endl;

    std::string seq = "ACGTACACTGGACTGAGTCT";
    std::string qual = "$$%%)*0)'%%%%)*0)'%%";
    btllib::PhredNtHash dna_phred_nthash(seq, 3, 20, 0, qual);

    std::string rna_seq = "ACGUACACUGGACUGAGUCU";
    btllib::PhredNtHash rna_phred_nthash(rna_seq, 3, 20, 0, qual);

    dna_phred_nthash.roll();
    rna_phred_nthash.roll();
    size_t i;
    for (i = 0; i < dna_phred_nthash.get_hash_num(); ++i) {
      TEST_ASSERT_EQ(dna_phred_nthash.hashes()[i], rna_phred_nthash.hashes()[i]);
    }
    TEST_ASSERT_EQ(i, 3);
  }

  return 0;
}

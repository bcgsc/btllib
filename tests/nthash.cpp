#include "btllib/nthash.hpp"
#include "btllib/seq.hpp"
#include "helpers.hpp"

#include <iostream>
#include <queue>
#include <stack>
#include <string>

int
main()
{

  {
    PRINT_TEST_NAME("k-mer hash values")

    std::string seq = "ACATGCATGCA";
    const unsigned k = 5;
    const unsigned h = 3;

    const std::vector<std::array<uint64_t, h>> hashes = {
      { 0x38cc00f940aebdae, 0xab7e1b110e086fc6, 0x11a1818bcfdd553 },
      { 0x603a48c5a11c794a, 0xe66016e61816b9c4, 0xc5b13cb146996ffe }
    };

    btllib::NtHash nthash(seq, h, k);
    nthash.roll();
    btllib::BlindNtHash blind(seq.data(), h, k);

    for (const auto& h_vals : hashes) {
      nthash.roll();
      TEST_ASSERT_ARRAY_EQ(h_vals, nthash.hashes(), h);
      blind.roll(seq[blind.get_pos() + k]);
      TEST_ASSERT_ARRAY_EQ(h_vals, blind.hashes(), h);
    }
  }

  {
    PRINT_TEST_NAME("k-mer rolling")

    std::string seq = "AGTCAGTC";
    unsigned h = 3;
    unsigned k = 4;

    btllib::NtHash nthash(seq, h, k);
    std::vector<uint64_t*> hashes;

    while (nthash.roll()) {
      uint64_t* h_vals = new uint64_t[h];
      std::copy(nthash.hashes(), nthash.hashes() + h, h_vals);
      hashes.push_back(h_vals);
    }

    TEST_ASSERT_EQ(hashes.size(), seq.length() - k + 1);

    // check same hash value for identical k-mers (first and last)
    TEST_ASSERT_ARRAY_EQ(hashes.front(), hashes.back(), h);
  }

  {
    PRINT_TEST_NAME("k-mer rolling vs. base hash values")

    std::string seq = "ACGTACACTGGACTGAGTCT";

    btllib::NtHash nthash(seq, 3, seq.size() - 2);
    /* 18-mers of kmer*/
    std::string kmer1 = seq.substr(0, 18);
    std::string kmer2 = seq.substr(1, 18);
    std::string kmer3 = seq.substr(2, 18);

    btllib::NtHash nthash_vector[] = {
      btllib::NtHash(kmer1, nthash.get_hash_num(), kmer1.size()),
      btllib::NtHash(kmer2, nthash.get_hash_num(), kmer2.size()),
      btllib::NtHash(kmer3, nthash.get_hash_num(), kmer3.size())
    };

    size_t i;
    for (i = 0; nthash.roll() && nthash_vector[i].roll(); ++i) {
      for (size_t j = 0; j < nthash.get_hash_num(); ++j) {
        TEST_ASSERT_EQ(nthash.hashes()[j], nthash_vector[i].hashes()[j]);
      }
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    PRINT_TEST_NAME("canonical hashing")

    std::string seq_f = "ACGTACACTGGACTGAGTCT";
    std::string seq_r = "AGACTCAGTCCAGTGTACGT";
    unsigned h = 3;

    btllib::NtHash nthash_f(seq_f, h, seq_f.size());
    btllib::NtHash nthash_r(seq_r, h, seq_r.size());

    nthash_f.roll();
    nthash_r.roll();
    TEST_ASSERT_EQ(nthash_f.get_hash_num(), nthash_r.get_hash_num())
    TEST_ASSERT_ARRAY_EQ(nthash_f.hashes(), nthash_r.hashes(), h)
  }

  {
    PRINT_TEST_NAME("k-mer back rolling")

    std::string seq = "ACTAGCTG";
    unsigned h = 3;
    unsigned k = 5;

    btllib::NtHash nthash(seq, h, k);
    std::stack<uint64_t*> hashes;

    while (nthash.roll()) {
      uint64_t* h_vals = new uint64_t[h];
      std::copy(nthash.hashes(), nthash.hashes() + h, h_vals);
      hashes.push(h_vals);
    }

    TEST_ASSERT_EQ(hashes.size(), seq.length() - k + 1)

    do {
      TEST_ASSERT_ARRAY_EQ(nthash.hashes(), hashes.top(), h)
      hashes.pop();
    } while (nthash.roll_back());
  }

  {
    PRINT_TEST_NAME("k-mer peeking")

    std::string seq = "ACTGATCAG";
    unsigned h = 3;
    unsigned k = 6;

    btllib::NtHash nthash(seq, h, k);
    nthash.roll();

    size_t steps = 3;
    while (steps--) {
      nthash.peek();
      uint64_t* h_peek = new uint64_t[h];
      std::copy(nthash.hashes(), nthash.hashes() + h, h_peek);
      nthash.peek(seq[nthash.get_pos() + k]);
      TEST_ASSERT_ARRAY_EQ(nthash.hashes(), h_peek, h);
      nthash.roll();
      TEST_ASSERT_ARRAY_EQ(nthash.hashes(), h_peek, h);
    }
  }

  {
    PRINT_TEST_NAME("skipping Ns")

    std::string seq = "ACGTACACTGGACTGAGTCT";
    std::string seq_with_ns = seq;

    TEST_ASSERT_GE(seq_with_ns.size(), 10)
    seq_with_ns[seq_with_ns.size() / 2] = 'N';
    seq_with_ns[seq_with_ns.size() / 2 + 1] = 'N';
    unsigned k = (seq.size() - 2) / 2 - 1;
    btllib::NtHash nthash(seq_with_ns, 3, k);

    std::vector<uint64_t> positions;
    for (size_t i = 0; i < seq_with_ns.size() / 2 - k + 1; i++) {
      positions.push_back(i);
    }
    for (size_t i = seq_with_ns.size() / 2 + 2; i < seq_with_ns.size() - k + 1;
         i++) {
      positions.push_back(i);
    }

    size_t i = 0;
    while (nthash.roll()) {
      TEST_ASSERT_EQ(nthash.get_pos(), positions[i])
      i++;
    }
    TEST_ASSERT_EQ(positions.size(), i)
  }

  {
    PRINT_TEST_NAME("base substitution")

    std::string seq = "ACGTACACTGGACTGAGTCT";
    std::string sub = "ACGCGCACTGGACTGAGTCT";

    btllib::NtHash nthash(seq, 3, seq.size());
    btllib::NtHash nthash_subbed(sub, 3, sub.size());

    nthash.roll();
    nthash.sub({ 3, 4 }, { 'C', 'G' });
    nthash_subbed.roll();
    TEST_ASSERT_EQ(nthash.get_hash_num(), nthash_subbed.get_hash_num());
    size_t i;
    for (i = 0; i < nthash.get_hash_num(); ++i) {
      TEST_ASSERT_EQ(nthash.hashes()[i], nthash_subbed.hashes()[i]);
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    PRINT_TEST_NAME("RNA")
    unsigned h = 3, k = 20;

    std::string seq = "ACGTACACTGGACTGAGTCTACGG";
    btllib::NtHash dna_nthash(seq, h, k);

    std::string rna_seq = "ACGUACACUGGACUGAGUCUACGG";
    btllib::NtHash rna_nthash(rna_seq, h, k);

    bool can_roll = true;
    while (can_roll) {
      can_roll = dna_nthash.roll();
      can_roll &= rna_nthash.roll();
      TEST_ASSERT_ARRAY_EQ(dna_nthash.hashes(), rna_nthash.hashes(), h);
    }
  }

  {
    PRINT_TEST_NAME("spaced seed hash values")

    std::string seq = "ACATGCATGCA";
    std::vector<std::string> seeds = { "11100111" };
    const unsigned k = seeds[0].length();
    const unsigned h = 3;

    const std::vector<std::array<uint64_t, h>> hashes = {
      { 0x10be4904ad8de5d, 0x3e29e4f4c991628c, 0x3f35c984b13feb20 },
      { 0x8200a7aa3eaf17c8, 0x344198402f4c2a9c, 0xb6423fe62e69c40c },
      { 0x3ce8adcbeaa56532, 0x162e91a4dbedbf11, 0x53173f786a031f45 }
    };

    btllib::SeedNtHash nthash(seq, seeds, h, k);

    for (const auto& h_vals : hashes) {
      nthash.roll();
      TEST_ASSERT_ARRAY_EQ(h_vals, nthash.hashes(), h);
    }
  }

  {
    PRINT_TEST_NAME("spaced seeds")

    std::string seq = "ACGTACACTGGACTGAGTCT";
    std::vector<std::string> seeds = { "111110000000011111",
                                       "111111100001111111" };

    /* Point mutations of k-mer */
    std::string seqM1 = "ACGTACACTTGACTGAGTCT";
    std::string seqM2 = "ACGTACACTGTACTGAGTCT";
    std::string seqM3 = "ACGTACACTGCACTGAGTCT";

    unsigned k = seq.size() - 2;
    TEST_ASSERT_EQ(k, seeds[0].size());
    TEST_ASSERT_EQ(k, seeds[1].size());

    btllib::SeedNtHash seed_nthash(seq, seeds, 2, k);
    btllib::SeedNtHash seed_nthashM1(seqM1, seeds, 2, k);
    btllib::SeedNtHash seed_nthashM2(seqM2, seeds, 2, k);
    btllib::SeedNtHash seed_nthashM3(seqM3, seeds, 2, k);

    std::vector<std::vector<uint64_t>> hashes;

    TEST_ASSERT_EQ(seed_nthash.get_hash_num(), seeds.size() * 2);

    size_t steps = 0;
    for (; seed_nthash.roll(); steps++) {
      TEST_ASSERT_EQ(seed_nthashM1.roll(), true);
      TEST_ASSERT_EQ(seed_nthashM2.roll(), true);
      TEST_ASSERT_EQ(seed_nthashM3.roll(), true);

      const std::string seq_sub = seq.substr(steps, k);
      const std::string seqM1_sub = seqM1.substr(steps, k);
      const std::string seqM2_sub = seqM2.substr(steps, k);
      const std::string seqM3_sub = seqM3.substr(steps, k);
      btllib::SeedNtHash seed_nthash_base(seq_sub, seeds, 2, k);
      btllib::SeedNtHash seed_nthashM1_base(seqM1_sub, seeds, 2, k);
      btllib::SeedNtHash seed_nthashM2_base(seqM2_sub, seeds, 2, k);
      btllib::SeedNtHash seed_nthashM3_base(seqM3_sub, seeds, 2, k);

      TEST_ASSERT_EQ(seed_nthash_base.roll(), true);
      TEST_ASSERT_EQ(seed_nthashM1_base.roll(), true);
      TEST_ASSERT_EQ(seed_nthashM2_base.roll(), true);
      TEST_ASSERT_EQ(seed_nthashM3_base.roll(), true);

      hashes.push_back({});
      for (size_t i = 0; i < seed_nthash.get_hash_num(); i++) {
        const auto hval = seed_nthash.hashes()[i];
        TEST_ASSERT_EQ(hval, seed_nthashM1.hashes()[i]);
        TEST_ASSERT_EQ(hval, seed_nthashM2.hashes()[i]);
        TEST_ASSERT_EQ(hval, seed_nthashM3.hashes()[i]);
        TEST_ASSERT_EQ(hval, seed_nthashM1_base.hashes()[i]);
        TEST_ASSERT_EQ(hval, seed_nthashM2_base.hashes()[i]);
        TEST_ASSERT_EQ(hval, seed_nthashM3_base.hashes()[i]);
        hashes.back().push_back(hval);
      }

      if (seed_nthash.get_pos() > 0) {
        seed_nthash.peek_back();
        for (size_t i = 0; i < seed_nthash.get_hash_num(); i++) {
          TEST_ASSERT_EQ(seed_nthash.hashes()[i], hashes[hashes.size() - 2][i]);
        }
        seed_nthash.peek_back(seq[seed_nthash.get_pos() - 1]);
        for (size_t i = 0; i < seed_nthash.get_hash_num(); i++) {
          TEST_ASSERT_EQ(seed_nthash.hashes()[i], hashes[hashes.size() - 2][i]);
        }
      }
    }
    TEST_ASSERT_EQ(seed_nthashM1.roll(), false);
    TEST_ASSERT_EQ(seed_nthashM2.roll(), false);
    TEST_ASSERT_EQ(seed_nthashM3.roll(), false);
    TEST_ASSERT_EQ(steps, seq.size() - k + 1);
  }

  {
    PRINT_TEST_NAME("spaced seed back roll")

    std::string seq = "ACTAGCTG";
    std::string seed = "110011";
    unsigned h = 3;
    unsigned k = seed.length();

    btllib::SeedNtHash nthash(seq, { seed }, h, k);
    std::stack<uint64_t*> hashes;

    while (nthash.roll()) {
      uint64_t* h_vals = new uint64_t[h];
      std::copy(nthash.hashes(), nthash.hashes() + h, h_vals);
      hashes.push(h_vals);
    }

    TEST_ASSERT_EQ(hashes.size(), seq.length() - k + 1)

    do {
      TEST_ASSERT_ARRAY_EQ(nthash.hashes(), hashes.top(), h)
      hashes.pop();
    } while (nthash.roll_back());
  }

  {
    PRINT_TEST_NAME("canonical hashing in spaced seeds")

    std::string seq_fwd = "CACTCGGCCACACACACACACACACACCCTCACACACACAAAACGCACAC";
    std::string seq_rev = "GTGTGCGTTTTGTGTGTGTGAGGGTGTGTGTGTGTGTGTGTGGCCGAGTG";

    std::vector<std::string> seeds = {
      "11011000001100101101011000011010110100110000011011",
      "01010000101001110100111011011100101110010100001010",
      "11100000100111010111000100100011101011100100000111",
      "01111000011000111101000011000010111100011000011110",
      "00111000011000111101000011000010111100011000011100",
      "00000000000000000000000011000000000000000000000000",
      "11111111111111111111111100111111111111111111111111",
      "11111111111111111111111111111111111111111111111111",
    };

    const unsigned h = 4;

    btllib::SeedNtHash h1(seq_fwd, seeds, h, seeds[0].length());
    btllib::SeedNtHash h2(seq_rev, seeds, h, seeds[0].length());

    bool can_roll = true;
    while (can_roll) {
      can_roll = h1.roll();
      can_roll &= h2.roll();
      TEST_ASSERT_ARRAY_EQ(h1.hashes(), h2.hashes(), h);
    }
  }

  {
    PRINT_TEST_NAME("copying SeedNtHash objects")

    std::string seq = "AACGTGACTACTGACTAGCTAGCTAGCTGATCGT";
    std::vector<std::string> seeds = { "111111111101111111111",
                                       "110111010010010111011" };

    const unsigned h = 4;

    btllib::SeedNtHash h1(seq, seeds, h, seeds[0].length());
    btllib::SeedNtHash h2(h1);

    bool can_roll = true;
    while (can_roll) {
      can_roll = h1.roll();
      can_roll &= h2.roll();
      TEST_ASSERT_ARRAY_EQ(h1.hashes(), h2.hashes(), h);
    }
  }

  {
    PRINT_TEST_NAME("BlindSeedNtHash")

    std::string seq = "ATGCTAGTAGCTGAC";
    std::vector<std::string> seeds = { "110011", "101101" };

    btllib::SeedNtHash h1(seq, seeds, 3, seeds[0].size());
    h1.roll();
    btllib::BlindSeedNtHash h2(seq.data(), seeds, 3, seeds[0].size());

    while (h1.roll()) {
      h2.roll(seq[h2.get_pos() + seeds[0].size()]);
      TEST_ASSERT_ARRAY_EQ(h1.hashes(), h2.hashes(), 6)
    }
  }

  {
    PRINT_TEST_NAME("BlindSeedNtHash roll back")

    std::string kmer = "ACCAGT";
    std::vector<std::string> seeds = { "110011", "101101" };

    btllib::BlindSeedNtHash h(kmer.data(), seeds, 3, seeds[0].size());
    h.roll('A');
    const auto hashes1 = h.hashes();
    h.roll_back('A');
    TEST_ASSERT_ARRAY_EQ(hashes1, h.hashes(), 3)
  }

  {
    PRINT_TEST_NAME("BlindSeedNtHash copy constructor")

    std::string seq = "ATGCTAGTAGCTGAC";
    std::vector<std::string> seeds = { "110011", "101101" };

    btllib::BlindSeedNtHash h1(seq.data(), seeds, 1, seeds[0].size());
    h1.roll('A');
    h1.roll('C');
    btllib::BlindSeedNtHash h2(h1);
    TEST_ASSERT_ARRAY_EQ(h1.hashes(), h2.hashes(), 2)
    h1.roll('G');
    h2.roll('G');
    TEST_ASSERT_ARRAY_EQ(h1.hashes(), h2.hashes(), 2)
    h1.roll('T');
    h2.roll('T');
    TEST_ASSERT_ARRAY_EQ(h1.hashes(), h2.hashes(), 2)
  }

  {
    PRINT_TEST_NAME("k-mer vs. full-care spaced seed hashing")
    const std::string seq = "ATGCTAGTAGCTGAC";
    const std::vector<std::string> seeds = { "11111" };
    const unsigned k = seeds[0].size();
    const unsigned h = 3;

    btllib::NtHash kmer(seq, h, k);
    btllib::SeedNtHash seed(seq, seeds, h, k);

    bool can_roll = true;
    while (can_roll) {
      can_roll = kmer.roll();
      can_roll |= seed.roll();
      TEST_ASSERT_ARRAY_EQ(kmer.hashes(), seed.hashes(), h)
    }
  }

  {
    PRINT_TEST_NAME("k-mer reset seq")
    const std::string seq1 = "ATGCTAGTAGCTGAC";
    const std::string seq2 = "GTGACTAGCTGACTG";
    const unsigned h = 3, k = 7;
    btllib::NtHash h1(seq1, h, k);
    h1.roll();
    h1.roll();
    h1.roll();
    h1.set_seq(seq2);
    btllib::NtHash h2(seq2, h, k);
    bool can_roll = true;
    while (can_roll) {
      can_roll = h1.roll();
      can_roll |= h2.roll();
      TEST_ASSERT_ARRAY_EQ(h1.hashes(), h2.hashes(), h)
    }
    TEST_ASSERT(!h1.roll())
    TEST_ASSERT(!h2.roll())
  }

  {
    PRINT_TEST_NAME("seed reset seq")
    const std::string seq1 = "ATGCTAGTAGCTGAC";
    const std::string seq2 = "GTGACTAGCTGACTG";
    const std::vector<std::string> seeds = { "11011" };
    const unsigned h = 3, k = seeds[0].size();
    btllib::SeedNtHash h1(seq1, seeds, h, k);
    h1.roll();
    h1.roll();
    h1.roll();
    h1.set_seq(seq2);
    btllib::SeedNtHash h2(seq2, seeds, h, k);
    bool can_roll = true;
    while (can_roll) {
      can_roll = h1.roll();
      can_roll |= h2.roll();
      TEST_ASSERT_ARRAY_EQ(h1.hashes(), h2.hashes(), h)
    }
    TEST_ASSERT(!h1.roll())
    TEST_ASSERT(!h2.roll())
  }

  return 0;
}

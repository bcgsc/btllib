#include "btllib/nthash.hpp"
#include "helpers.hpp"

#include <iostream>
#include <string>

int
main()
{
  const std::string seq = "ACGTACACTGGACTGAGTCT";

  {
    std::cerr << "Testing single kmer hash values" << std::endl;
    btllib::NtHash nthash(seq, 3, seq.size());

    const std::vector<uint64_t> hashes = { 10434435546371013747UL,
                                           16073887395445158014UL,
                                           8061578976118370557UL };

    nthash.roll();
    TEST_ASSERT_EQ(nthash.get_hash_num(), hashes.size());
    size_t i;
    for (i = 0; i < nthash.get_hash_num(); ++i) {
      TEST_ASSERT_EQ(nthash.hashes()[i], hashes[i]);
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    std::cerr << "Testing rolling" << std::endl;
    btllib::NtHash nthash(seq, 3, seq.size() - 2);

    const std::vector<uint64_t> hashes = {
      13163759238790597551UL, 12621194261766804115UL, 7338209520836137437UL,
      1791600248699339671UL,  11191929811376415288UL, 12983529916756035783UL,
      5425935920907136555UL,  5015085639900805921UL,  10441021516240327793UL
    };

    size_t steps = 0;
    while (nthash.roll()) {
      for (size_t i = 0; i < nthash.get_hash_num(); ++i) {
        TEST_ASSERT_EQ(nthash.hashes()[i],
                       hashes[nthash.get_pos() * nthash.get_hash_num() + i]);
      }
      if (nthash.get_pos() > 0) {
        nthash.peek_back();
        for (size_t i = 0; i < nthash.get_hash_num(); ++i) {
          TEST_ASSERT_EQ(
            nthash.hashes()[i],
            hashes[(nthash.get_pos() - 1) * nthash.get_hash_num() + i]);
        }
        nthash.peek_back(seq[nthash.get_pos() - 1]);
        for (size_t i = 0; i < nthash.get_hash_num(); ++i) {
          TEST_ASSERT_EQ(
            nthash.hashes()[i],
            hashes[(nthash.get_pos() - 1) * nthash.get_hash_num() + i]);
        }
      }
      steps++;
    }
    TEST_ASSERT_EQ(steps, nthash.get_pos() + 1)
    TEST_ASSERT_EQ(steps, 3)
  }

  {
    std::cerr << "Testing rolling backwards" << std::endl;
    unsigned k = seq.size() - 2;
    btllib::NtHash nthash(seq, 3, k, seq.size() - k);

    const std::vector<uint64_t> hashes = {
      13163759238790597551UL, 12621194261766804115UL, 7338209520836137437UL,
      1791600248699339671UL,  11191929811376415288UL, 12983529916756035783UL,
      5425935920907136555UL,  5015085639900805921UL,  10441021516240327793UL
    };

    size_t steps = 3;
    while (nthash.roll_back()) {
      for (size_t i = 0; i < nthash.get_hash_num(); ++i) {
        TEST_ASSERT_EQ(nthash.hashes()[i],
                       hashes[nthash.get_pos() * nthash.get_hash_num() + i]);
      }
      if (nthash.get_pos() < seq.size() - k) {
        nthash.peek();
        for (size_t i = 0; i < nthash.get_hash_num(); ++i) {
          TEST_ASSERT_EQ(
            nthash.hashes()[i],
            hashes[(nthash.get_pos() + 1) * nthash.get_hash_num() + i]);
        }
        nthash.peek(seq[nthash.get_pos() + k]);
        for (size_t i = 0; i < nthash.get_hash_num(); ++i) {
          TEST_ASSERT_EQ(
            nthash.hashes()[i],
            hashes[(nthash.get_pos() + 1) * nthash.get_hash_num() + i]);
        }
      }
      steps--;
    }
    TEST_ASSERT_EQ(nthash.get_pos(), 0)
    TEST_ASSERT_EQ(steps, 0)
  }

  {
    std::cerr << "Testing rolling backwards with spaced seeds" << std::endl;
    unsigned k = seq.size() - 2;
    std::vector<std::string> seeds = { "111110000000011111",
                                       "111111100001111111" };
    btllib::SeedNtHash seed_nthash(seq, seeds, 2, k, seq.size() - k);

    const std::vector<uint64_t> hashes = {
      817502191096265638ULL,   7002589659100769832ULL,  1859879627538729626ULL,
      1503737960973537115ULL,  11201337262633184115ULL, 8091243882347916556ULL,
      13566405847740282550ULL, 9923089871377977073ULL,  1451958618186721217ULL,
      3905378141411987905ULL,  6006262340350391097ULL,  8505903736152692418ULL
    };

    size_t steps = 3;
    while (seed_nthash.roll_back()) {
      for (size_t i = 0; i < seed_nthash.get_hash_num(); ++i) {
        TEST_ASSERT_EQ(
          seed_nthash.hashes()[i],
          hashes[seed_nthash.get_pos() * seed_nthash.get_hash_num() + i]);
      }
      if (seed_nthash.get_pos() < seq.size() - k) {
        seed_nthash.peek();
        for (size_t i = 0; i < seed_nthash.get_hash_num(); ++i) {
          TEST_ASSERT_EQ(
            seed_nthash.hashes()[i],
            hashes[(seed_nthash.get_pos() + 1) * seed_nthash.get_hash_num() +
                   i]);
        }
        seed_nthash.peek(seq[seed_nthash.get_pos() + k]);
        for (size_t i = 0; i < seed_nthash.get_hash_num(); ++i) {
          TEST_ASSERT_EQ(
            seed_nthash.hashes()[i],
            hashes[(seed_nthash.get_pos() + 1) * seed_nthash.get_hash_num() +
                   i]);
        }
      }
      steps--;
    }
    TEST_ASSERT_EQ(seed_nthash.get_pos(), 0)
    TEST_ASSERT_EQ(steps, 0)
  }

  {
    std::cerr << "Testing skipping Ns" << std::endl;
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
    std::cerr << "Testing base substitution" << std::endl;
    btllib::NtHash nthash(seq, 3, seq.size());
    std::string seq_subbed = "ACGCGCACTGGACTGAGTCT";
    btllib::NtHash nthash_subbed(seq_subbed, 3, seq_subbed.size());

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
    std::cerr << "Testing reverse complement" << std::endl;
    /* Reverse complement of kmer*/
    std::string rc_seq = "AGACTCAGTCCAGTGTACGT";

    btllib::NtHash nthash(seq, 3, seq.size());
    btllib::NtHash nthash_rc(rc_seq, 3, rc_seq.size());

    nthash.roll();
    nthash_rc.roll();
    TEST_ASSERT_EQ(nthash.get_hash_num(), nthash_rc.get_hash_num());
    size_t i;
    for (i = 0; i < nthash.get_hash_num(); ++i) {
      TEST_ASSERT_EQ(nthash.hashes()[i], nthash_rc.hashes()[i]);
    }
    TEST_ASSERT_EQ(i, 3);
  }

  {
    std::cerr << "Testing rolling vs ntbase hash values" << std::endl;
    btllib::NtHash nthash(seq, 3, seq.size() - 2);

    /* 18-mers of kmer*/
    std::string kmer1 = "ACGTACACTGGACTGAGT";
    std::string kmer2 = "CGTACACTGGACTGAGTC";
    std::string kmer3 = "GTACACTGGACTGAGTCT";

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
    std::cerr << "Testing spaced seeds" << std::endl;
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
    std::cerr << "Testing RNA" << std::endl;
    btllib::NtHash dna_nthash(seq, 3, 20);

    std::string rna_seq = "ACGUACACUGGACUGAGUCU";
    btllib::NtHash rna_nthash(rna_seq, 3, 20);

    dna_nthash.roll();
    rna_nthash.roll();
    size_t i;
    for (i = 0; i < dna_nthash.get_hash_num(); ++i) {
      TEST_ASSERT_EQ(dna_nthash.hashes()[i], rna_nthash.hashes()[i]);
    }
    TEST_ASSERT_EQ(i, 3);
  }

  return 0;
}
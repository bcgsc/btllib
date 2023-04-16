#include "btllib/aahash.hpp"

#include "helpers.hpp"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

int main()
{
  const std::string seq = "RITMLYTIRITMLYTI";
  const unsigned k = 8, h = 3;
  PRINT_TEST_NAME("Level 1 hashing")
  btllib::AAHash aahash(seq, h, k, 1);
  unsigned num_kmers = seq.size() - k + 2;
  aahash.roll();
  --num_kmers;
  uint64_t first_hash = aahash.hashes()[0];
  while (--num_kmers)
  {
    TEST_ASSERT(aahash.roll());
  }
  TEST_ASSERT(!aahash.roll());
  uint64_t last_hash = aahash.hashes()[0];
  TEST_ASSERT_EQ(first_hash, last_hash)

  PRINT_TEST_NAME("Level 2 hashing")
  btllib::AAHash aahash2(seq, h, k, 2);
  num_kmers = seq.size() - k + 2;
  aahash2.roll();
  --num_kmers;
  first_hash = aahash2.hashes()[0];
  while (--num_kmers)
  {
    TEST_ASSERT(aahash2.roll());
  }
  TEST_ASSERT(!aahash2.roll());
  last_hash = aahash2.hashes()[0];
  TEST_ASSERT_EQ(first_hash, last_hash)

  PRINT_TEST_NAME("Level 3 hashing")
  btllib::AAHash aahash3(seq, h, k, 3);
  num_kmers = seq.size() - k + 2;
  aahash3.roll();
  --num_kmers;
  first_hash = aahash3.hashes()[0];
  while (--num_kmers)
  {
    TEST_ASSERT(aahash3.roll());
  }
  TEST_ASSERT(!aahash3.roll());
  last_hash = aahash3.hashes()[0];
  TEST_ASSERT_EQ(first_hash, last_hash)

  PRINT_TEST_NAME("Test seed parse")

  std::string string_seed1 = "0123";
  std::string string_seed2 = "1230";
  std::string string_seed3 = "2301";
  std::string string_seed4 = "3012";

  std::vector<std::string> seeds = {string_seed1, string_seed2, string_seed3, string_seed4};

  btllib::SpacedSeed unsigned_seed1 = {0, 1, 2, 3};
  btllib::SpacedSeed unsigned_seed2 = {1, 2, 3, 0};
  btllib::SpacedSeed unsigned_seed3 = {2, 3, 0, 1};
  btllib::SpacedSeed unsigned_seed4 = {3, 0, 1, 2};

  std::vector<btllib::SpacedSeed> seeds_ss = {unsigned_seed1, unsigned_seed2, unsigned_seed3, unsigned_seed4};

  std::vector<btllib::SpacedSeed> parsed_seeds = btllib::aa_parse_seeds(seeds);

  TEST_ASSERT_EQ(parsed_seeds.size(), seeds_ss.size());

  for (size_t i = 0; i < parsed_seeds.size(); ++i)
  {
    for (size_t j = 0; j < parsed_seeds[i].size(); ++j)
    {
      TEST_ASSERT_EQ(parsed_seeds[i][j], seeds_ss[i][j]);
    }
  }

  PRINT_TEST_NAME("multi-level hashing")
  std::string string_seed5 = "01233210";
  std::string string_seed6 = "12300321";
  std::string string_seed7 = "23011032";
  std::string string_seed8 = "30122103";
  std::vector<std::string> seeds2 = {string_seed5, string_seed6, string_seed7, string_seed8};
  std::vector<btllib::SpacedSeed> parsed_seeds2 = btllib::aa_parse_seeds(seeds2);
  btllib::SeedAAHash seedaahash(seq, parsed_seeds2, h, k);
  btllib::SeedAAHash seedaahash2(seq, parsed_seeds2, h, k);
  num_kmers = seq.size() - k + 2;
  seedaahash.roll();
  seedaahash2.roll();
  --num_kmers;
  while (--num_kmers)
  {
    TEST_ASSERT(seedaahash.roll());
  }
  TEST_ASSERT(!seedaahash.roll());
  for (size_t i = 0; i < seeds2.size() * h; ++i)
  {
    TEST_ASSERT_EQ(seedaahash.hashes()[i], seedaahash2.hashes()[i]);
  }
  TEST_ASSERT_EQ(first_hash, last_hash)

  PRINT_TEST_NAME("Nonequivalence")
  std::string seq2 = "ALGASQMYWMCDYDPY";
  btllib::AAHash aahash4(seq2, h, k, 1);
  btllib::AAHash aahash5(seq, h, k, 1);
  num_kmers = seq.size() - k + 1;
  aahash4.roll();
  aahash5.roll();
  TEST_ASSERT_NE(aahash4.hashes()[0], aahash5.hashes()[0]);
  --num_kmers;
  while (--num_kmers)
  {
    TEST_ASSERT_NE(aahash4.hashes()[0], aahash5.hashes()[0]);
  }

  PRINT_TEST_NAME("Level 2 equivalence")
  std::string seqT = "T";
  std::string seqS = "S";
  btllib::AAHash aahashlvl2T(seqT, 1, 1, 2);
  btllib::AAHash aahashlvl2S(seqS, 1, 1, 2);
  aahashlvl2T.roll();
  aahashlvl2S.roll();
  TEST_ASSERT_EQ(aahashlvl2T.hashes()[0], aahashlvl2S.hashes()[0]);
  std::string seqD = "D";
  std::string seqE = "E";
  btllib::AAHash aahashlvl2D(seqD, 1, 1, 2);
  btllib::AAHash aahashlvl2E(seqE, 1, 1, 2);
  aahashlvl2D.roll();
  aahashlvl2E.roll();
  TEST_ASSERT_EQ(aahashlvl2D.hashes()[0], aahashlvl2E.hashes()[0]);
  std::string seqQ = "Q";
  std::string seqK = "K";
  std::string seqR = "R";
  btllib::AAHash aahashlvl2Q(seqQ, 1, 1, 2);
  btllib::AAHash aahashlvl2K(seqK, 1, 1, 2);
  btllib::AAHash aahashlvl2R(seqR, 1, 1, 2);
  aahashlvl2Q.roll();
  aahashlvl2K.roll();
  aahashlvl2R.roll();
  TEST_ASSERT_EQ(aahashlvl2Q.hashes()[0], aahashlvl2K.hashes()[0]);
  TEST_ASSERT_EQ(aahashlvl2Q.hashes()[0], aahashlvl2R.hashes()[0]);
  std::string seqV = "V";
  std::string seqI = "I";
  std::string seqL = "L";
  std::string seqM = "M";
  btllib::AAHash aahashlvl2V(seqV, 1, 1, 2);
  btllib::AAHash aahashlvl2I(seqI, 1, 1, 2);
  btllib::AAHash aahashlvl2L(seqL, 1, 1, 2);
  btllib::AAHash aahashlvl2M(seqM, 1, 1, 2);
  aahashlvl2V.roll();
  aahashlvl2I.roll();
  aahashlvl2L.roll();
  aahashlvl2M.roll();
  TEST_ASSERT_EQ(aahashlvl2V.hashes()[0], aahashlvl2I.hashes()[0]);
  TEST_ASSERT_EQ(aahashlvl2V.hashes()[0], aahashlvl2L.hashes()[0]);
  TEST_ASSERT_EQ(aahashlvl2V.hashes()[0], aahashlvl2M.hashes()[0]);
  std::string seqW = "W";
  std::string seqF = "F";
  std::string seqY = "Y";
  btllib::AAHash aahashlvl2W(seqW, 1, 1, 2);
  btllib::AAHash aahashlvl2F(seqF, 1, 1, 2);
  btllib::AAHash aahashlvl2Y(seqY, 1, 1, 2);
  aahashlvl2W.roll();
  aahashlvl2F.roll();
  aahashlvl2Y.roll();
  TEST_ASSERT_EQ(aahashlvl2W.hashes()[0], aahashlvl2F.hashes()[0]);
  TEST_ASSERT_EQ(aahashlvl2W.hashes()[0], aahashlvl2Y.hashes()[0]);

  std::string seqCGATNDQVWHP = "CGATNDQVWHP";
  std::string seqCGASNEKIFHP = "CGASNEKIFHP";
  btllib::AAHash aahashlvl2CGATNDQVWHP(seqCGATNDQVWHP, 1, 11, 2);
  btllib::AAHash aahashlvl2CGASNEKIFHP(seqCGASNEKIFHP, 1, 11, 2);
  aahashlvl2CGATNDQVWHP.roll();
  aahashlvl2CGASNEKIFHP.roll();
  TEST_ASSERT_EQ(aahashlvl2CGATNDQVWHP.hashes()[0], aahashlvl2CGASNEKIFHP.hashes()[0]);

  PRINT_TEST_NAME("Level 3 equivalence")
  std::string seqA = "A";
  btllib::AAHash aahashlvl3T(seqT, 1, 1, 3);
  btllib::AAHash aahashlvl3S(seqS, 1, 1, 3);
  btllib::AAHash aahashlvl3A(seqA, 1, 1, 3);
  aahashlvl3T.roll();
  aahashlvl3S.roll();
  aahashlvl3A.roll();
  TEST_ASSERT_EQ(aahashlvl3T.hashes()[0], aahashlvl3S.hashes()[0]);
  TEST_ASSERT_EQ(aahashlvl3T.hashes()[0], aahashlvl3A.hashes()[0]);
  std::string seqN = "N";
  btllib::AAHash aahashlvl3D(seqD, 1, 1, 3);
  btllib::AAHash aahashlvl3E(seqE, 1, 1, 3);
  btllib::AAHash aahashlvl3N(seqN, 1, 1, 3);
  aahashlvl3D.roll();
  aahashlvl3E.roll();
  aahashlvl3N.roll();
  TEST_ASSERT_EQ(aahashlvl3D.hashes()[0], aahashlvl3E.hashes()[0]);
  TEST_ASSERT_EQ(aahashlvl3D.hashes()[0], aahashlvl3N.hashes()[0]);
  btllib::AAHash aahashlvl3Q(seqQ, 1, 1, 3);
  btllib::AAHash aahashlvl3K(seqK, 1, 1, 3);
  btllib::AAHash aahashlvl3R(seqR, 1, 1, 3);
  aahashlvl3Q.roll();
  aahashlvl3K.roll();
  aahashlvl3R.roll();
  TEST_ASSERT_EQ(aahashlvl3Q.hashes()[0], aahashlvl3K.hashes()[0]);
  TEST_ASSERT_EQ(aahashlvl3Q.hashes()[0], aahashlvl3R.hashes()[0]);
  btllib::AAHash aahashlvl3V(seqV, 1, 1, 3);
  btllib::AAHash aahashlvl3I(seqI, 1, 1, 3);
  btllib::AAHash aahashlvl3L(seqL, 1, 1, 3);
  btllib::AAHash aahashlvl3M(seqM, 1, 1, 3);
  aahashlvl3V.roll();
  aahashlvl3I.roll();
  aahashlvl3L.roll();
  aahashlvl3M.roll();
  TEST_ASSERT_EQ(aahashlvl3V.hashes()[0], aahashlvl3I.hashes()[0]);
  TEST_ASSERT_EQ(aahashlvl3V.hashes()[0], aahashlvl3L.hashes()[0]);
  TEST_ASSERT_EQ(aahashlvl3V.hashes()[0], aahashlvl3M.hashes()[0]);
  btllib::AAHash aahashlvl3W(seqW, 1, 1, 3);
  btllib::AAHash aahashlvl3F(seqF, 1, 1, 3);
  btllib::AAHash aahashlvl3Y(seqY, 1, 1, 3);
  aahashlvl3W.roll();
  aahashlvl3F.roll();
  aahashlvl3Y.roll();
  TEST_ASSERT_EQ(aahashlvl3W.hashes()[0], aahashlvl3F.hashes()[0]);
  TEST_ASSERT_EQ(aahashlvl3W.hashes()[0], aahashlvl3Y.hashes()[0]);

  std::string seqCGANQVWHP = "CGANQVWHP";
  std::string seqCGTDKIFHP = "CGTDKIFHP";
  btllib::AAHash aahashlvl3CGANQVWHP(seqCGANQVWHP, 1, 9, 3);
  btllib::AAHash aahashlvl3CGTDKIFHP(seqCGTDKIFHP, 1, 9, 3);
  aahashlvl3CGANQVWHP.roll();
  aahashlvl3CGTDKIFHP.roll();
  TEST_ASSERT_EQ(aahashlvl3CGANQVWHP.hashes()[0], aahashlvl3CGTDKIFHP.hashes()[0]);

  return 0;
}
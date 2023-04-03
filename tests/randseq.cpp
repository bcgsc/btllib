#include "btllib/randseq.hpp"
#include "btllib/status.hpp"

#include <algorithm>
#include <cctype>
#include <unordered_set>

using namespace btllib;

void test_DNASequence() {
  RandSeq generator(RandSeq::SeqType::DNA);
  std::string seq = generator.generate(100);

  check_error(seq.size() != 100, "DNA sequence length is incorrect.");
  for (char c : seq) {
    check_error(!(c == 'A' || c == 'C' || c == 'G' || c == 'T'), "Invalid character in DNA sequence.");
  }
  log_info("DNA sequence test passed.");
}

void test_RNASequence() {
  RandSeq generator(RandSeq::SeqType::RNA);
  std::string seq = generator.generate(100);

  check_error(seq.size() != 100, "RNA sequence length is incorrect.");
  for (char c : seq) {
    check_error(!(c == 'A' || c == 'C' || c == 'G' || c == 'U'), "Invalid character in RNA sequence.");
  }
  log_info("RNA sequence test passed.");
}

void test_PROTEINSequence() {
  RandSeq generator(RandSeq::SeqType::PROTEIN);
  std::string seq = generator.generate(100);

  check_error(seq.size() != 100, "Protein sequence length is incorrect.");
  for (char c : seq) {
    check_error(!std::isupper(c) ||
    !(std::unordered_set<char>{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}.count(c)),
               "Invalid character in Protein sequence.");
  }
  log_info("Protein sequence test passed.");
}

void test_Masking() {
  RandSeq generator(RandSeq::SeqType::DNA, RandSeq::Masking::SOFT);
  std::string seq = generator.generate(100);

  check_error(seq.size() != 100, "DNA sequence with soft masking length is incorrect.");
  for (char c : seq) {
    check_error(!(std::toupper(c) == 'A' || std::toupper(c) == 'C' || std::toupper(c) == 'G' || std::toupper(c) == 'T'),
               "Invalid character in DNA sequence with soft masking.");
  }
  log_info("Soft masking test passed.");
}

int main() {
  test_DNASequence();
  test_RNASequence();
  test_PROTEINSequence();
  test_Masking();
  return 0;
}

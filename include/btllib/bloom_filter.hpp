#ifndef BTLLIB_BLOOM_FILTER_HPP
#define BTLLIB_BLOOM_FILTER_HPP

#include "status.hpp"

#include "vendor/cpptoml.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace btllib {

class BloomFilter
{

public:
  BloomFilter() {}
  BloomFilter(size_t size, unsigned hash_num);
  BloomFilter(const std::string& path);

  ~BloomFilter();

  void insert(const std::vector<uint64_t>& hashes);
  void insert(const uint64_t* hashes);

  bool contains(const std::vector<uint64_t>& hashes);
  bool contains(const uint64_t* hashes);

  size_t get_size();
  uint64_t get_pop_cnt();
  unsigned get_hash_num();
  double get_fpr();

  void write(const std::string& path);

private:
  static unsigned pop_cnt_byte(unsigned char x);

  static constexpr unsigned BITS_IN_BYTE = 8; // NOLINT
  static constexpr unsigned char BIT_MASKS[] = {
    // NOLINT
    0x01, 0x02, 0x04, 0x08, // NOLINT
    0x10, 0x20, 0x40, 0x80  // NOLINT
  };
  inline static const char* MAGIC_HEADER = "BTLBloomFilter_v1";

  unsigned char* bitarray = nullptr;
  size_t size = 0; // In bytes
  unsigned hash_num = 0;
};

inline BloomFilter::BloomFilter(size_t size, unsigned hash_num)
  : size(std::ceil(size / sizeof(long)) * sizeof(long))
  , hash_num(hash_num)
{
  delete[] bitarray;
  bitarray = new unsigned char[this->size];
#pragma omp parallel for
  for (size_t i = 0; i < size / sizeof(long); i++) {
    *(((long*)bitarray) + i) = 0;
  }
}

inline BloomFilter::BloomFilter(const std::string& path)
{
  std::ifstream file(path);

  std::string magic_with_brackets = std::string("[") + MAGIC_HEADER + "]";

  std::string line;
  std::getline(file, line);
  if (line != magic_with_brackets) {
    log_error(
      std::string("Magic string does not match (likely version mismatch)\n") +
      "Your magic string:        " + line + "\n" +
      "BloomFilter magic string: " + magic_with_brackets);
    std::exit(EXIT_FAILURE);
  }

  /* Read bloom filter line by line until it sees "[HeaderEnd]"
  which is used to mark the end of the header section and
  assigns the header to a char array*/
  std::string toml_buffer(line + '\n');
  bool header_end_found = false;
  while (std::getline(file, line)) {
    toml_buffer.append(line + '\n');
    if (line == "[HeaderEnd]") {
      header_end_found = true;
      break;
    }
  }
  if (!header_end_found) {
    log_error("Pre-built bloom filter does not have the correct header end.");
    std::exit(EXIT_FAILURE);
  }

  // Send the char array to a stringstream for the cpptoml parser to parse
  std::istringstream toml_stream(toml_buffer);
  cpptoml::parser toml_parser(toml_stream);
  auto header_config = toml_parser.parse();

  // Obtain header values from toml parser and assign them to class members
  auto table = header_config->get_table(MAGIC_HEADER);
  size = *table->get_as<size_t>("size");
  hash_num = *table->get_as<unsigned>("hash_num");

  delete[] bitarray;
  bitarray = new unsigned char[size];
  file.read((char*)bitarray, size);
}

inline BloomFilter::~BloomFilter()
{
  delete[] bitarray;
}

inline void
BloomFilter::insert(const std::vector<uint64_t>& hashes)
{
  insert(hashes.data());
}

inline void
BloomFilter::insert(const uint64_t* hashes)
{
  for (unsigned i = 0; i < hash_num; ++i) {
    auto normalized = hashes[i] % size;
    __sync_or_and_fetch(&(bitarray[normalized / BITS_IN_BYTE]),
                        BIT_MASKS[normalized % BITS_IN_BYTE]);
  }
}

inline bool
BloomFilter::contains(const std::vector<uint64_t>& hashes)
{
  return contains(hashes.data());
}

inline bool
BloomFilter::contains(const uint64_t* hashes)
{
  for (unsigned i = 0; i < hash_num; ++i) {
    auto normalized = hashes[i] % size;
    auto mask = BIT_MASKS[normalized % BITS_IN_BYTE];
    if (!bool(bitarray[normalized / BITS_IN_BYTE] & mask)) {
      return false;
    }
  }
  return true;
}

inline size_t
BloomFilter::get_size()
{
  return size;
}

inline uint64_t
BloomFilter::get_pop_cnt()
{
  uint64_t pop_cnt = 0;
#pragma omp parallel for reduction(+ : pop_cnt)
  for (size_t i = 0; i < size; ++i) {
    pop_cnt += pop_cnt_byte(bitarray[i]);
  }
  return pop_cnt;
}

inline unsigned
BloomFilter::get_hash_num()
{
  return hash_num;
}

inline double
BloomFilter::get_fpr()
{
  return std::pow(double(get_pop_cnt()) / double(size), double(hash_num));
}

inline void
BloomFilter::write(const std::string& path)
{
  std::ofstream file(path.c_str(), std::ios::out | std::ios::binary);

  /* Initialize cpptoml root table
    Note: Tables and fields are unordered
    Ordering of table is maintained by directing the table
    to the output stream immediately after completion  */
  auto root = cpptoml::make_table();

  /* Initialize bloom filter section and insert fields
      and output to ostream */
  auto header = cpptoml::make_table();
  header->insert("size", size);
  header->insert("hash_num", hash_num);
  root->insert(MAGIC_HEADER, header);
  file << *root << "[HeaderEnd]\n";

  file.write((char*)bitarray, size);
}

inline unsigned
BloomFilter::pop_cnt_byte(unsigned char x)
{
  return ((0x876543210 >>                                              // NOLINT
           (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >> // NOLINT
          ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2)) &  // NOLINT
         0xf;                                                          // NOLINT
}

} // namespace btllib

#endif
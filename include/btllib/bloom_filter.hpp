#ifndef BTLLIB_BLOOM_FILTER_HPP
#define BTLLIB_BLOOM_FILTER_HPP

#include "rolling_hash.hpp"
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
  BloomFilter(size_t bytes, unsigned hash_num);
  BloomFilter(const std::string& path);

  void insert(const std::vector<uint64_t>& hashes);
  void insert(const uint64_t* hashes);

  bool contains(const std::vector<uint64_t>& hashes);
  bool contains(const uint64_t* hashes);

  size_t get_bytes() const { return bytes; }
  uint64_t get_pop_cnt() const;
  unsigned get_hash_num() const { return hash_num; }
  double get_fpr() const;

  void write(const std::string& path);

private:
  std::vector<unsigned char> bytearray;
  size_t bytes = 0;
  unsigned hash_num = 0;
};

class CountingBloomFilter
{

public:
  CountingBloomFilter() {}
  CountingBloomFilter(size_t bytes, unsigned hash_num);
  CountingBloomFilter(const std::string& path);

  void insert(const std::vector<uint64_t>& hashes);
  void insert(const uint64_t* hashes);

  unsigned contains(const std::vector<uint64_t>& hashes);
  unsigned contains(const uint64_t* hashes);

  size_t get_bytes() const { return bytes; }
  uint64_t get_pop_cnt() const;
  unsigned get_hash_num() const { return hash_num; }
  double get_fpr() const;

  void write(const std::string& path);

private:
  std::vector<unsigned char> bytearray;
  size_t bytes = 0;
  unsigned hash_num = 0;
};

class KmerBloomFilter
{

public:
  KmerBloomFilter(unsigned k, size_t bytes, unsigned hash_num = 4);

  void insert(const std::string& seq);
  void insert(const char* seq, size_t seq_len);

  unsigned contains(const std::string& seq);
  unsigned contains(const char* seq, size_t seq_len);

private:
  unsigned k;
  BloomFilter bf;
};

class KmerCountingBloomFilter
{

public:
  KmerCountingBloomFilter(unsigned k, size_t bytes, unsigned hash_num = 4);

  void insert(const std::string& seq);
  void insert(const char* seq, size_t seq_len);

  unsigned contains(const std::string& seq);
  unsigned contains(const char* seq, size_t seq_len);

private:
  unsigned k;
  CountingBloomFilter bf;
};

static const unsigned BITS_IN_BYTE = 8; // NOLINT
static const unsigned char BIT_MASKS[BITS_IN_BYTE] = {
  // NOLINT
  0x01, 0x02, 0x04, 0x08, // NOLINT
  0x10, 0x20, 0x40, 0x80  // NOLINT
};
static const char* BLOOM_FILTER_MAGIC_HEADER = "BTLBloomFilter_v2";
static const char* COUNTING_BLOOM_FILTER_MAGIC_HEADER =
  "BTLCountingBloomFilter_v2";

inline unsigned
pop_cnt_byte(unsigned char x)
{
  return ((0x876543210 >>                                              // NOLINT
           (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >> // NOLINT
          ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2)) &  // NOLINT
         0xf;                                                          // NOLINT
}

inline BloomFilter::BloomFilter(size_t bytes, unsigned hash_num)
  : bytes(std::ceil(bytes / sizeof(long)) * sizeof(long))
  , hash_num(hash_num)
{
  bytearray.resize(bytes);
}

inline BloomFilter::BloomFilter(const std::string& path)
{
  std::ifstream file(path);

  std::string magic_with_brackets =
    std::string("[") + BLOOM_FILTER_MAGIC_HEADER + "]";

  std::string line;
  std::getline(file, line);
  if (line != magic_with_brackets) {
    log_error(
      std::string("Magic string does not match (likely version mismatch)\n") +
      "Your magic string:\t" + line + "\n" + "BloomFilter magic string:\t" +
      magic_with_brackets);
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
  auto table = header_config->get_table(BLOOM_FILTER_MAGIC_HEADER);
  bytes = *table->get_as<size_t>("bytes");
  hash_num = *table->get_as<unsigned>("hash_num");

  bytearray.resize(bytes);
  file.read((char*)bytearray.data(), bytes);
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
    auto normalized = hashes[i] % bytes;
    __sync_or_and_fetch(&(bytearray[normalized / BITS_IN_BYTE]),
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
    auto normalized = hashes[i] % bytes;
    auto mask = BIT_MASKS[normalized % BITS_IN_BYTE];
    if (!bool(bytearray[normalized / BITS_IN_BYTE] & mask)) {
      return false;
    }
  }
  return true;
}

inline uint64_t
BloomFilter::get_pop_cnt() const
{
  uint64_t pop_cnt = 0;
#pragma omp parallel for reduction(+ : pop_cnt)
  for (size_t i = 0; i < bytes; ++i) {
    pop_cnt += pop_cnt_byte(bytearray[i]);
  }
  return pop_cnt;
}

inline double
BloomFilter::get_fpr() const
{
  return std::pow(double(get_pop_cnt()) / double(bytes), double(hash_num));
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
  header->insert("bytes", bytes);
  header->insert("hash_num", hash_num);
  root->insert(BLOOM_FILTER_MAGIC_HEADER, header);
  file << *root << "[HeaderEnd]\n";

  file.write((char*)bytearray.data(), bytes);
}

inline CountingBloomFilter::CountingBloomFilter(size_t bytes, unsigned hash_num)
  : bytes(bytes)
  , hash_num(hash_num)
{}

inline CountingBloomFilter::CountingBloomFilter(const std::string& path) {}

inline void
CountingBloomFilter::insert(const std::vector<uint64_t>& hashes)
{
  insert(hashes.data());
}

inline void
CountingBloomFilter::insert(const uint64_t* hashes)
{
  (void)hashes;
}

inline unsigned
CountingBloomFilter::contains(const std::vector<uint64_t>& hashes)
{
  return contains(hashes.data());
}

inline unsigned
CountingBloomFilter::contains(const uint64_t* hashes)
{
  (void)hashes;
  return 0;
}

inline uint64_t
CountingBloomFilter::get_pop_cnt() const
{
  return 0;
}

inline double
CountingBloomFilter::get_fpr() const
{
  return 0;
}

inline void
CountingBloomFilter::write(const std::string& path)
{}

inline KmerBloomFilter::KmerBloomFilter(unsigned k,
                                        size_t bytes,
                                        unsigned hash_num)
  : k(k)
  , bf(bytes, hash_num)
{}

inline void
KmerBloomFilter::insert(const std::string& seq)
{
  insert(seq.c_str(), seq.size());
}

inline void
KmerBloomFilter::insert(const char* seq, size_t seq_len)
{
  RollingHash rolling_hash(seq, seq_len, k, bf.get_hash_num());
  while (rolling_hash.roll()) {
    bf.insert(rolling_hash.hashes());
  }
}

inline unsigned
KmerBloomFilter::contains(const std::string& seq)
{
  return contains(seq.c_str(), seq.size());
}
inline unsigned
KmerBloomFilter::contains(const char* seq, size_t seq_len)
{
  unsigned count = 0;
  RollingHash rolling_hash(seq, seq_len, k, bf.get_hash_num());
  while (rolling_hash.roll()) {
    if (bf.contains(rolling_hash.hashes())) {
      count++;
    }
  }
  return count;
}

inline KmerCountingBloomFilter::KmerCountingBloomFilter(unsigned k,
                                                        size_t bytes,
                                                        unsigned hash_num)
  : k(k)
  , bf(bytes, hash_num)
{}

inline void
KmerCountingBloomFilter::insert(const std::string& seq)
{
  insert(seq.c_str(), seq.size());
}

inline void
KmerCountingBloomFilter::insert(const char* seq, size_t seq_len)
{
  RollingHash rolling_hash(seq, seq_len, k, bf.get_hash_num());
  while (rolling_hash.roll()) {
    bf.insert(rolling_hash.hashes());
  }
}

inline unsigned
KmerCountingBloomFilter::contains(const std::string& seq)
{
  return contains(seq.c_str(), seq.size());
}

inline unsigned
KmerCountingBloomFilter::contains(const char* seq, size_t seq_len)
{
  unsigned count = 0;
  RollingHash rolling_hash(seq, seq_len, k, bf.get_hash_num());
  while (rolling_hash.roll()) {
    count += bf.contains(rolling_hash.hashes());
  }
  return count;
}

} // namespace btllib

#endif
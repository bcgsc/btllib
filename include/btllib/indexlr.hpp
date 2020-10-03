#ifndef BTLLIB_INDEXLR_HPP
#define BTLLIB_INDEXLR_HPP

#include "nthash.hpp"
#include "order_queue.hpp"
#include "seq_reader.hpp"
#include "status.hpp"
#include "util.hpp"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace btllib {

static const unsigned BUFFER_SIZE = 512;
static const unsigned BLOCK_SIZE = 64;

static const std::string END_TOKEN = "01a44a7c-d9b6-4994-a49d-b1811f92c6e1";

// TODO: Allow multiple Indexlr objects to be instantiated (by assigning ID to
// each instance / indexing static members based on ID)
class Indexlr
{

public:
  enum Flag
  {
    ID = 0,
    NO_ID = 1,
    BX = 2,
    NO_BX = 0,
    SEQ = 4,
    NO_SEQ = 0
  };

  bool output_id() const { return bool(~flags & NO_ID); }
  bool output_bx() const { return bool(flags & BX); }
  bool output_seq() const { return bool(flags & SEQ); }

  struct Read
  {
    Read() {}

    Read(size_t num, std::string id, std::string comment, std::string seq)
      : num(num)
      , id(std::move(id))
      , comment(std::move(comment))
      , seq(std::move(seq))
    {}

    size_t num = 0;
    std::string id;
    std::string comment;
    std::string seq;
  };

  struct Minimizer
  {
    uint64_t hash1, hash2;
    size_t pos;
    std::string seq;
  };

  using HashedKmer = Minimizer;

  struct Record
  {
    Record() {}

    Record(size_t num,
           std::string id,
           std::string barcode,
           std::vector<Minimizer> minimizers)
      : num(num)
      , id(std::move(id))
      , barcode(std::move(barcode))
      , minimizers(std::move(minimizers))
    {}

    size_t num = 0;
    std::string id;
    std::string barcode;
    std::vector<Minimizer> minimizers;

    operator bool() const { return !id.empty() || !barcode.empty(); }
  };

  Record get_minimizers();

  Indexlr(const std::string& seqfile,
          size_t k,
          size_t w,
          unsigned flags = 0,
          unsigned thread = 1);

private:
  static std::string extract_barcode(std::string comment);
  std::vector<HashedKmer> hash_kmers(const std::string& seq, size_t k) const;
  static std::vector<Minimizer> minimize_hashed_kmers(
    const std::vector<HashedKmer>& hashed_kmers,
    size_t w);

  const std::string& seqfile;
  const size_t k, w;
  unsigned flags;
  unsigned threads;

  std::atomic<bool> fasta{ false };
  OrderQueueSPMC<Read, BUFFER_SIZE, BLOCK_SIZE> input_queue;
  OrderQueueMPSC<Record, BUFFER_SIZE, BLOCK_SIZE> output_queue;

  class Worker
  {
  public:
    void start() { t = std::thread(do_work, this); }
    void join() { t.join(); }

    virtual ~Worker() {}

    Worker& operator=(const Worker& worker) = delete;
    Worker& operator=(Worker&& worker) = delete;

  protected:
    Worker(Indexlr& indexlr)
      : indexlr(indexlr)
    {}

    Worker(const Worker& worker)
      : Worker(worker.indexlr)
    {}
    Worker(Worker&& worker) noexcept
      : Worker(worker.indexlr)
    {}

    Indexlr& indexlr;

    virtual void work() = 0;
    static void do_work(Worker* worker) { worker->work(); }

    std::thread t;
  };

  class InputWorker : public Worker
  {
  public:
    InputWorker(Indexlr& indexlr)
      : Worker(indexlr)
    {}

    InputWorker(const InputWorker& worker)
      : InputWorker(worker.indexlr)
    {}
    InputWorker(InputWorker&& worker) noexcept
      : InputWorker(worker.indexlr)
    {}

    InputWorker& operator=(const InputWorker& worker) = delete;
    InputWorker& operator=(InputWorker&& worker) = delete;

    void work() override;
  };

  class MinimizeWorker : public Worker
  {
  public:
    MinimizeWorker(Indexlr& indexlr)
      : Worker(indexlr)
    {}

    MinimizeWorker(const MinimizeWorker& worker)
      : MinimizeWorker(worker.indexlr)
    {}
    MinimizeWorker(MinimizeWorker&& worker) noexcept
      : MinimizeWorker(worker.indexlr)
    {}

    MinimizeWorker& operator=(const MinimizeWorker& worker) = delete;
    MinimizeWorker& operator=(MinimizeWorker&& worker) = delete;

    void work() override;
  };
};

inline Indexlr::Indexlr(const std::string& seqfile,
                        const size_t k,
                        const size_t w,
                        const unsigned flags,
                        const unsigned threads)
  : seqfile(seqfile)
  , k(k)
  , w(w)
  , flags(flags)
  , threads(threads)
{
  InputWorker iw(*this);
  iw.start();
  auto mw = std::vector<MinimizeWorker>(threads, MinimizeWorker(*this));
  for (auto& worker : mw) {
    worker.start();
  }
  for (auto& worker : mw) {
    worker.join();
  }
  iw.join();
}

// Minimerize a sequence: Find the minimizers of a vector of hash values
// representing a sequence.
/* Algorithm
v is a vector of non-negative integers
w is the window size
Invariants
    0 <  w <= v.size() - 1
    0 <= l <= r <= v.size() - 1
Initial conditions
    M    = NIL       Final set of minimizers, empty initially
    min  = -1        Minimum element
    i    = -1        Index of minimum element
    prev = -1        Index of previous minimum element
    l    = 0         Index of left end of window
    r    = l + w - 1 Index of right end of window
Computation
At each window, if the previous minimum is out of scope, find the new,
right-most, minimum or else, check with only the right-most element to determine
if that is the new minimum. A minimizer is added to the final vector only if
it's index has changed. for each window of v bounded by [l, r] if (i < l) i =
index of minimum element in [l, r], furthest from l. else if (v[r] <= v[i]) i =
r min = v[i] if (i != prev) { prev = i M <- M + m
    }
    l = l + 1        Move window's left bound by one element
    r = l + w - 1    Set window's right bound
}*/

inline std::string
Indexlr::extract_barcode(std::string comment)
{
  const static std::string BARCODE_PREFIX = "BX:Z:";
  if (starts_with(comment, BARCODE_PREFIX)) {
    auto pos = comment.find(' ');
    if (pos != std::string::npos) {
      comment.erase(pos);
    }
    comment.erase(0, BARCODE_PREFIX.size());
  } else {
    comment = "NA";
  }
  return std::move(comment);
}

inline std::vector<Indexlr::HashedKmer>
Indexlr::hash_kmers(const std::string& seq, const size_t k) const
{
  std::vector<HashedKmer> hashed_kmers;
  if (seq.size() < k) {
    return {};
  }
  hashed_kmers.reserve(seq.size() - k + 1);
  for (NtHash nh(seq, k, 2); nh.roll();) {
    hashed_kmers.push_back(
      Minimizer({ nh.hashes()[0],
                  nh.hashes()[1],
                  nh.get_pos(),
                  output_seq() ? seq.substr(nh.get_pos(), k) : "" }));
  }
  return hashed_kmers;
}

inline std::vector<Indexlr::Minimizer>
Indexlr::minimize_hashed_kmers(
  const std::vector<Indexlr::HashedKmer>& hashed_kmers,
  const size_t w)
{
  if (hashed_kmers.size() < w) {
    return {};
  }
  std::vector<Minimizer> minimizers;
  minimizers.reserve(2 * hashed_kmers.size() / w);
  int i = -1, prev = -1;
  auto first_it = hashed_kmers.begin();
  auto min_it = hashed_kmers.end();
  for (auto left_it = first_it; left_it < hashed_kmers.end() - w + 1;
       ++left_it) {
    auto right_it = left_it + w;
    if (i < left_it - first_it) {
      // Use of operator '<=' returns the minimum that is furthest from left.
      min_it = std::min_element(
        left_it, right_it, [](const HashedKmer& a, const HashedKmer& b) {
          return a.hash1 <= b.hash1;
        });
    } else if (right_it[-1].hash1 <= min_it->hash1) {
      min_it = right_it - 1;
    }
    i = min_it - first_it;
    if (i > prev) {
      prev = i;
      minimizers.push_back(*min_it);
    }
  }
  return minimizers;
}

inline Indexlr::Record
Indexlr::get_minimizers()
{
  thread_local static decltype(output_queue)::Block block;
  if (block.count <= block.current) {
    output_queue.read(block);
    if (block.count < block.current) {
      block = decltype(output_queue)::Block();
      return Record();
    }
  }
  auto& record = block.data[block.current];
  if (record.id == END_TOKEN) {
    block = decltype(output_queue)::Block();
    return Record();
  }
  block.current++;
  return std::move(record);
}

inline void
Indexlr::InputWorker::work()
{
  SeqReader reader(indexlr.seqfile);
  if (reader.get_format() == reader.FASTA) {
    indexlr.fasta = true;
  } else {
    indexlr.fasta = false;
  }

  decltype(indexlr.input_queue)::Block block;
  size_t current_block_num = 0;
  SeqReader::Record record;
  Read read;
  while ((record = reader.read())) {
    block.data[block.count++] = Read(record.num,
                                     std::move(record.name),
                                     std::move(record.comment),
                                     std::move(record.seq));
    if (block.count == BLOCK_SIZE) {
      block.num = current_block_num++;
      indexlr.input_queue.write(block);
      block.count = 0;
    }
  }
  for (unsigned i = 0; i < indexlr.threads; i++) {
    block.data[block.count++] = Read(0, END_TOKEN, "", "");
    if (block.count == BLOCK_SIZE || i == indexlr.threads - 1) {
      block.num = current_block_num++;
      indexlr.input_queue.write(block);
      block.current = 0;
      block.count = 0;
    }
  }
}

inline void
Indexlr::MinimizeWorker::work()
{
  decltype(indexlr.input_queue)::Block input_block;
  decltype(indexlr.output_queue)::Block output_block;

  for (;;) {
    if (input_block.current == input_block.count) {
      indexlr.input_queue.read(input_block);
    }
    Read& read = input_block.data[input_block.current++];
    if (read.id == END_TOKEN) {
      output_block.data[output_block.count++] =
        Record(read.num, std::move(read.id), std::move(read.comment), {});
      output_block.num = input_block.num;
      indexlr.output_queue.write(output_block);
      break;
    }
    Record record;
    record.num = read.num;
    if (indexlr.output_id()) {
      record.id = std::move(read.id);
    }
    if (indexlr.output_bx()) {
      record.barcode = indexlr.extract_barcode(read.comment);
    }

    check_warning(read.seq.size() < indexlr.k,
                  "Indexlr: skipped seq " + std::to_string(read.num) +
                    " on line " +
                    std::to_string(read.num * (indexlr.fasta ? 2 : 4) + 2) +
                    "; k (" + std::to_string(indexlr.k) + ") > seq length (" +
                    std::to_string(read.seq.size()) + ")");

    decltype(indexlr.hash_kmers(read.seq, indexlr.k)) hashed_kmers;
    if (read.seq.size() >= indexlr.k) {
      hashed_kmers = indexlr.hash_kmers(read.seq, indexlr.k);

      check_warning(
        indexlr.w > hashed_kmers.size(),
        "Indexlr: skipped seq " + std::to_string(read.num) + " on line " +
          std::to_string(read.num * (indexlr.fasta ? 2 : 4) + 2) + "; w (" +
          std::to_string(indexlr.w) + ") > # of hashes (" +
          std::to_string(hashed_kmers.size()) + ")");
    }

    if (indexlr.w <= hashed_kmers.size() && read.seq.size() >= indexlr.k) {
      record.minimizers =
        indexlr.minimize_hashed_kmers(hashed_kmers, indexlr.w);
    } else {
      record.minimizers = {};
    }

    output_block.data[output_block.count++] = std::move(record);
    if (output_block.count == BLOCK_SIZE) {
      output_block.num = input_block.num;
      indexlr.output_queue.write(output_block);
      output_block.current = 0;
      output_block.count = 0;
    }
  }
}

} // namespace btllib

#endif
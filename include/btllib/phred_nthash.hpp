#ifndef BTLLIB_PHRED_NTHASH_HPP
#define BTLLIB_PHRED_NTHASH_HPP

#include "btllib/nthash.hpp"
#include "btllib/util.hpp"

#include <vector>
#include <string>

namespace btllib
{
    /**
     * NtHash with Phred score filtering.
     */
    class PhredNtHash : private NtHash
    {
    public:
        /**
         * Constructor for PhredNtHash.
         * @param seq Sequence to hash.
         * @param seq_len Length of `seq`.
         * @param hash_num Number of hashes to compute.
         * @param k Length of k-mer.
         * @param phred_min Minimum Phred score for a base to be included in the hash.
         * @param quality_string String of Phred scores for each base in `seq`.
         * @param pos Position to start hashing from.
         * @return PhredNtHash object.
         */
        PhredNtHash(const char *seq,
                    size_t seq_len,
                    unsigned hash_num,
                    unsigned k,
                    size_t phred_min, const std::string &quality_string, size_t pos = 0);

        /**
         * Constructor for PhredNtHash.
         * @param seq Sequence to hash.
         * @param hash_num Number of hashes to compute.
         * @param k Length of k-mer.
         * @param phred_min Minimum Phred score for a base to be included in the hash.
         * @param quality_string String of Phred scores for each base in `seq`.
         * @param pos Position to start hashing from.
         * @return PhredNtHash object.
         */
        PhredNtHash(const std::string &seq, unsigned hash_num, unsigned k,
                    size_t phred_min, const std::string &quality_string, size_t pos = 0);

        /**
         * Constructor for PhredNtHash.
         * @param seq Sequence to hash.
         * @param seq_len Length of `seq`.
         * @param hash_num Number of hashes to compute.
         * @param k Length of k-mer.
         * @param phred_min Minimum Phred score for a base to be included in the hash.
         * @param quality_string String of Phred scores for each base in `seq`.
         * @param pos Position to start hashing from.
         * @return PhredNtHash object.
         */
        PhredNtHash(const char *seq,
                    size_t seq_len,
                    unsigned hash_num,
                    unsigned k,
                    size_t phred_min, const char *quality_string, size_t pos = 0);

        /**
         * Constructor for PhredNtHash.
         * @param seq Sequence to hash.
         * @param hash_num Number of hashes to compute.
         * @param k Length of k-mer.
         * @param phred_min Minimum Phred score for a base to be included in the hash.
         * @param quality_string String of Phred scores for each base in `seq`.
         * @param pos Position to start hashing from.
         * @return PhredNtHash object.
         */
        PhredNtHash(const std::string &seq, unsigned hash_num, unsigned k,
                    size_t phred_min, const char *quality_string, size_t pos = 0);

        /**
         * Roll the hash forward by one base. If the next k-mer contains a base with
         * a Phred score below `phred_min`, the hash will be rolled forward until a
         * k-mer with all bases with a Phred score above `phred_min` is found.
         * @return True if successful, false if the end of the sequence is reached or
         * no k-mer with all bases with a Phred score above `phred_min` is found.
         */
        bool roll();
        /**
         * Roll the hash backward by one base. If the previous k-mer contains a base
         * with a Phred score below `phred_min`, the hash will be rolled backward
         * until a k-mer with all bases with a Phred score above `phred_min` is found.
         * @return True if successful, false if the start of the sequence is reached or
         * no k-mer with all bases with a Phred score above `phred_min` is found.
         */
        bool roll_back();
        const uint64_t *hashes() const { return NtHash::hashes(); }

        /**
         * Get the position of last hashed k-mer or the k-mer to be hashed if roll()
         * has never been called on this NtHash object.
         */
        size_t get_pos() const { return NtHash::get_pos(); }
        bool forward() const { return NtHash::forward(); }
        unsigned get_hash_num() const { return NtHash::get_hash_num(); }
        unsigned get_k() const { return NtHash::get_k(); }
        size_t get_seq_len() const { return NtHash::get_seq_len(); }
        uint64_t get_forward_hash() const { return NtHash::get_forward_hash(); }
        uint64_t get_reverse_hash() const { return NtHash::get_reverse_hash(); }

    private:
        const char *qual_seq;
        size_t phred_min;
        RangeMinimumQuery<std::string> rmq;
    };

} // namespace btllib

#endif // BTLLIB_PHRED_NTHASH_HPP
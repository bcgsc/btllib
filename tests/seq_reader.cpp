#include "../include/btl/seq_reader.hpp"

#include <iostream>

int main() {
    std::string seq, qual;

    const char* seqs[] = { "ACTG", "CATG" };
    const char* quals[] = { "!@^&", "(#&$" };

    int i;

    // Test FASTA
    btl::SeqReader reader_fasta("../tests/input.fa");
    assert(reader_fasta.get_format() == btl::SeqReader::Format::FASTA);

    i = 0;
	while (reader_fasta >> seq) {
        assert(seq == seqs[i]);
        assert(qual.empty());

        i++;
	}

    // Test FASTQ
    btl::SeqReader reader_fastq("../tests/input.fq");
    assert(reader_fastq.get_format() == btl::SeqReader::Format::FASTQ);

    i = 0;
	while (reader_fastq >> seq) {
        qual = reader_fastq.get_qual();

        assert(seq == seqs[i]);
        assert(qual == quals[i]);

        i++;
	}

    // Test SAM
    btl::SeqReader reader_sam("../tests/input.sam");
    assert(reader_sam.get_format() == btl::SeqReader::Format::SAM);

    i = 0;
	while (reader_sam >> seq) {
        qual = reader_sam.get_qual();

        assert(seq == seqs[i]);
        assert(qual == quals[i]);

        i++;
	}

	return 0;
}
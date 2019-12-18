#include "../include/btl/seq_reader.hpp"

#include <iostream>

int main() {
    std::string seq, qual;

    const char* seqs[] = { "ACTG", "TGCA" };
    const char* quals[] = { "!@^&", "(#&$" };

    int i;

    // Test FASTA
    btl::SeqReader reader_fasta("../tests/input.fa");
    assert(reader_fasta.get_format() == btl::SeqReader::Format::FASTA);

    i = 0;
	while (reader_fasta.read()) {
        seq = reader_fasta.seq();

        assert(seq == seqs[i]);
        assert(qual.empty());

        i++;
	}
    assert(i == 2);

    // Test FASTQ
    btl::SeqReader reader_fastq("../tests/input.fq");
    assert(reader_fastq.get_format() == btl::SeqReader::Format::FASTQ);

    i = 0;
	while (reader_fastq.read()) {
        seq = reader_fastq.seq();
        qual = reader_fastq.qual();

        assert(seq == seqs[i]);
        assert(qual == quals[i]);

        i++;
	}
    assert(i == 2);

    // Test SAM
    btl::SeqReader reader_sam("../tests/input.sam");
    assert(reader_sam.get_format() == btl::SeqReader::Format::SAM);

    i = 0;
	while (reader_sam.read()) {
        seq = reader_sam.seq();
        qual = reader_sam.qual();

        assert(seq == seqs[i]);
        assert(qual == quals[i]);

        i++;
	}
    assert(i == 2);

    // Test GFA2
    btl::SeqReader reader_gfa2("../tests/input.gfa2");
    assert(reader_gfa2.get_format() == btl::SeqReader::Format::GFA2);

    i = 0;
	while (reader_gfa2.read()) {
        seq = reader_gfa2.seq();
        qual = reader_gfa2.qual();

        assert(seq == seqs[i]);
        assert(qual.empty());

        i++;
	}
    assert(i == 2);

	return 0;
}
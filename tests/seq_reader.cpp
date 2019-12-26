#include "../include/btllib/seq_reader.hpp"

#include <cassert>
#include <string>

int main() {
    std::string seq, qual;

    const char* seqs[] = { "ACTG", "TGCA" };
    const char* quals[] = { "!@^&", "(#&$" };

    int i;

    // Test FASTA
    btllib::SeqReader reader_fasta("../tests/input.fa.gz");
    assert(reader_fasta.get_format() == btllib::SeqReader::Format::FASTA);

    i = 0;
	while (reader_fasta.read()) {
        seq = reader_fasta.seq();

        assert(seq == seqs[i]);
        assert(qual.empty());

        i++;
	}
    assert(i == 2);

    // Test FASTQ
    btllib::SeqReader reader_fastq("../tests/input.fq.tar.xz");
    assert(reader_fastq.get_format() == btllib::SeqReader::Format::FASTQ);

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
    btllib::SeqReader reader_sam("../tests/input.sam");
    assert(reader_sam.get_format() == btllib::SeqReader::Format::SAM);

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
    btllib::SeqReader reader_gfa2("../tests/input.gfa2");
    assert(reader_gfa2.get_format() == btllib::SeqReader::Format::GFA2);

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
#include "../include/btllib/seq_reader.hpp"

#include <cassert>
#include <string>

int main() {
    const char* seqs[] = { "ACTG", "TGCA" };
    const char* quals[] = { "!@^&", "(#&$" };

    int i;
    btllib::SeqReader::Record record;

    // Test FASTA
    btllib::SeqReader reader_fasta("../tests/input.fa.gz");
    assert(reader_fasta.get_format() == btllib::SeqReader::Format::FASTA);

    i = 0;
	while (record = reader_fasta.read()) {
        assert(record.seq == seqs[i]);
        assert(record.qual.empty());

        i++;
	}
    assert(i == 2);

    // Test FASTQ
    btllib::SeqReader reader_fastq("../tests/input.fq.tar.xz");
    assert(reader_fastq.get_format() == btllib::SeqReader::Format::FASTQ);

    i = 0;
	while (record = reader_fastq.read()) {
        assert(record.seq == seqs[i]);
        assert(record.qual == quals[i]);

        i++;
	}
    assert(i == 2);

    // Test SAM
    btllib::SeqReader reader_sam("../tests/input.bam");
    assert(reader_sam.get_format() == btllib::SeqReader::Format::SAM);

    i = 0;
	while (record = reader_sam.read()) {
        assert(record.seq == seqs[i]);
        assert(record.qual == quals[i]);

        i++;
	}
    assert(i == 2);

    // Test GFA2
    btllib::SeqReader reader_gfa2("../tests/input.gfa2");
    assert(reader_gfa2.get_format() == btllib::SeqReader::Format::GFA2);

    i = 0;
	while (record = reader_gfa2.read()) {
        assert(record.seq == seqs[i]);
        assert(record.qual.empty());

        i++;
	}
    assert(i == 2);

	return 0;
}
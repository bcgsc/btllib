#include "../include/btllib/seq_reader.hpp"
#include "../include/btllib/seq_writer.hpp"

#include <chrono>
#include <thread>

int main() {
    const char* names[] = { "1", "2" };
    const char* comments[] = { "comment1", "comment2" };
    const char* seqs[] = { "ACTG", "TGCA" };
    const char* quals[] = { "!@^&", "(#&$" };

    // Test FASTA
    btllib::SeqWriter writer_fasta("output.fa", btllib::SeqWriter::FASTA);
    for (int i = 0; i < 2; i++) {
        writer_fasta.write(names[i], comments[i], seqs[i], "");
    }
    writer_fasta.close();
    
    btllib::SeqReader reader_fasta("output.fa");
    assert(reader_fasta.get_format() == btllib::SeqReader::Format::FASTA);

    int j;
    btllib::SeqReader::Record record;

    j = 0;
	while (record = reader_fasta.read()) {
        assert(record.name == names[j]);
        assert(record.comment == comments[j]);
        assert(record.seq == seqs[j]);
        assert(record.qual.empty());

        j++;
	}
    assert(j == 2);

    // Test FASTQ
    btllib::SeqWriter writer_fastq("output.fq.bz2", btllib::SeqWriter::FASTQ);
    for (int i = 0; i < 2; i++) {
        writer_fastq.write(names[i], comments[i], seqs[i], quals[i]);
    }
    writer_fastq.close();
    
    btllib::SeqReader reader_fastq("output.fq.bz2");
    assert(reader_fastq.get_format() == btllib::SeqReader::Format::FASTQ);

    j = 0;
	while (record = reader_fastq.read()) {
        assert(record.name == names[j]);
        assert(record.comment == comments[j]);
        assert(record.seq == seqs[j]);
        assert(record.qual == quals[j]);

        j++;
	}
    assert(j == 2);

	return 0;
}
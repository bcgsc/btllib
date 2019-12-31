#include "../include/btllib/seq_reader.hpp"
#include "../include/btllib/seq_writer.hpp"

#include <chrono>
#include <thread>

int main() {
    std::string name, comment, seq, qual;

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

    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    
    btllib::SeqReader reader_fasta("output.fa");
    assert(reader_fasta.get_format() == btllib::SeqReader::Format::FASTA);

    int j;

    j = 0;
	while (reader_fasta.read()) {
        name = reader_fasta.name();
        comment = reader_fasta.comment();
        seq = reader_fasta.seq();

        assert(name == names[j]);
        assert(comment == comments[j]);
        assert(seq == seqs[j]);
        assert(qual.empty());

        j++;
	}
    assert(j == 2);

    // Test FASTQ
    btllib::SeqWriter writer_fastq("output.fq.bz2", btllib::SeqWriter::FASTQ);
    for (int i = 0; i < 2; i++) {
        writer_fastq.write(names[i], comments[i], seqs[i], quals[i]);
    }
    writer_fastq.close();

    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    
    btllib::SeqReader reader_fastq("output.fq.bz2");
    assert(reader_fastq.get_format() == btllib::SeqReader::Format::FASTQ);

    j = 0;
	while (reader_fastq.read()) {
        name = reader_fastq.name();
        comment = reader_fastq.comment();
        seq = reader_fastq.seq();
        qual = reader_fastq.qual();

        assert(name == names[j]);
        assert(comment == comments[j]);
        assert(seq == seqs[j]);
        assert(qual == quals[j]);

        j++;
	}
    assert(j == 2);

    // Meson deletes test files, so wait a bit for the pipes to close
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

	return 0;
}
#include "../include/btllib/seq_reader.hpp"
#include "../include/btllib/seq_writer.hpp"

#include "helpers.hpp"

#include <fstream>
#include <cstdio>

int main() {
    const char* names[] = { "1", "2" };
    const char* comments[] = { "comment1", "comment2" };
    const char* seqs[] = { "ACTG", "TGCA" };
    const char* quals[] = { "!@^&", "(#&$" };
    std::string filename;

    for (int iteration = 0; iteration < 3; iteration++) {
        std::cerr << "Iteration " << iteration + 1 << std::endl;

        // Test FASTA
        filename = get_random_name(64);
        std::cerr << "Test FASTA" << std::endl;
        btllib::SeqWriter writer_fasta(filename, btllib::SeqWriter::FASTA);
        for (int i = 0; i < 2; i++) {
            writer_fasta.write(names[i], comments[i], seqs[i], "");
        }
        writer_fasta.close();
        
        btllib::SeqReader reader_fasta(filename);
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

        reader_fasta.close();
        std::remove(filename.c_str());

        // Test FASTQ
        filename = get_random_name(64) + ".bz2";
        std::cerr << "Test FASTQ" << std::endl;
        btllib::SeqWriter writer_fastq(filename, btllib::SeqWriter::FASTQ);
        for (int i = 0; i < 2; i++) {
            writer_fastq.write(names[i], comments[i], seqs[i], quals[i]);
        }
        writer_fastq.close();
        
        btllib::SeqReader reader_fastq(filename);
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

        reader_fastq.close();
        std::remove(filename.c_str());

        // Test larger randomly generated file
        std::cerr << "Test random file" << std::endl;
        std::vector<std::string> generated_names;
        std::vector<std::string> generated_comments;
        std::vector<std::string> generated_seqs;
        std::vector<std::string> generated_quals;
        filename = get_random_name(64) + ".gz.xz";
        btllib::SeqWriter random_seqs(filename, btllib::SeqWriter::FASTQ, false);
        for (int s = 0; s < 500; s++) {
            std::string name, comment, seq, qual;

            name = get_random_name(10);
            comment = get_random_name(20);
            seq = get_random_sequence(200);
            qual = get_random_name(200);

            random_seqs.write(name, comment, seq, qual);

            generated_names.push_back(name);
            generated_comments.push_back(comment);
            generated_seqs.push_back(seq);
            generated_quals.push_back(qual);
        }
        random_seqs.close();

        btllib::SeqReader random_reader(filename);
        for (j = 0; record = random_reader.read(); j++) {
            assert(record.name == generated_names[j]);
            assert(record.comment == generated_comments[j]);
            assert(record.seq == generated_seqs[j]);
            assert(record.qual == generated_quals[j]);
        }
        assert(j == 500);

        random_reader.close();
        std::remove(filename.c_str());
    }

	return 0;
}
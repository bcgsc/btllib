#include "../include/btllib/seq_reader.hpp"

#include "helpers.hpp"

#include <cassert>
#include <string>
#include <fstream>
#include <cstdio>

int main() {
    const char* seqs[] = { "ACTG", "TGCA" };
    const char* quals[] = { "!@^&", "(#&$" };

    for (int iteration = 0; iteration < 3; iteration++) {
        std::cerr << "Iteration " << iteration + 1 << std::endl;

        int i;
        btllib::SeqReader::Record record;

        // Test FASTA
        std::cerr << "Test FASTA" << std::endl;
        btllib::SeqReader reader_fasta("../tests/input.fa.gz.bz2");
        assert(reader_fasta.get_format() == btllib::SeqReader::Format::FASTA);

        i = 0;
        while (record = reader_fasta.read()) {
            assert(record.seq == seqs[i]);
            assert(record.qual.empty());

            i++;
        }
        assert(i == 2);

        // Test FASTQ
        std::cerr << "Test FASTQ" << std::endl;
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
        std::cerr << "Test SAM" << std::endl;
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
        std::cerr << "Test GFA2" << std::endl;
        btllib::SeqReader reader_gfa2("../tests/input.gfa2");
        assert(reader_gfa2.get_format() == btllib::SeqReader::Format::GFA2);

        i = 0;
        while (record = reader_gfa2.read()) {
            assert(record.seq == seqs[i]);
            assert(record.qual.empty());

            i++;
        }
        assert(i == 2);

        // Test larger randomly generated file
        std::cerr << "Test random file" << std::endl;
        std::vector<std::string> generated_names;
        std::vector<std::string> generated_comments;
        std::vector<std::string> generated_seqs;
        std::vector<std::string> generated_quals;
        auto random_filename = get_random_name(64);
        std::ofstream random_seqs(random_filename);
        for (int s = 0; s < 500; s++) {
            std::string name, comment, seq, qual;

            name = get_random_name(10);
            comment = get_random_name(20);
            seq = get_random_sequence(200);
            qual = get_random_name(200);

            random_seqs << '@' << name << ' ' << comment << '\n'
              << seq << "\n+\n" << qual << '\n';

            generated_names.push_back(name);
            generated_comments.push_back(comment);
            generated_seqs.push_back(seq);
            generated_quals.push_back(qual);
        }
        random_seqs.close();

        btllib::SeqReader random_reader(random_filename);
        for (i = 0; record = random_reader.read(); i++) {
            assert(record.name == generated_names[i]);
            assert(record.comment == generated_comments[i]);
            assert(record.seq == generated_seqs[i]);
            assert(record.qual == generated_quals[i]);
        }
        assert(i == 500);
        
        random_reader.close();
        std::remove(random_filename.c_str());
    }

	return 0;
}
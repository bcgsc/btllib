#include "../include/btllib/seq_reader.hpp"
#include "../include/btllib/seq_writer.hpp"

#include <chrono>
#include <thread>
#include <fstream>
#include <random>

int main() {
    const char* names[] = { "1", "2" };
    const char* comments[] = { "comment1", "comment2" };
    const char* seqs[] = { "ACTG", "TGCA" };
    const char* quals[] = { "!@^&", "(#&$" };

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count() + 9999999);
    std::uniform_int_distribution<int> distribution_actg(0, 3);
    std::uniform_int_distribution<int> distribution_alphabet(65, 90);
    auto gen_random_actg = std::bind(distribution_actg, generator);
    auto gen_random_alphabet = std::bind(distribution_alphabet, generator);
    std::string random_filename;

    for (int iteration = 0; iteration < 3; iteration++) {
        std::cerr << "Iteration " << iteration + 1 << std::endl;

        random_filename.clear();
        for (int n = 0; n < 64; n++) {
            random_filename += char(gen_random_alphabet());
        }
        // Test FASTA
        std::cerr << "Test FASTA" << std::endl;
        btllib::SeqWriter writer_fasta(random_filename, btllib::SeqWriter::FASTA);
        for (int i = 0; i < 2; i++) {
            writer_fasta.write(names[i], comments[i], seqs[i], "");
        }
        writer_fasta.close();
        
        btllib::SeqReader reader_fasta(random_filename);
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

        random_filename.clear();
        for (int n = 0; n < 64; n++) {
            random_filename += char(gen_random_alphabet());
        }
        // Test FASTQ
        std::cerr << "Test FASTQ" << std::endl;
        btllib::SeqWriter writer_fastq(random_filename + ".bz2", btllib::SeqWriter::FASTQ);
        for (int i = 0; i < 2; i++) {
            writer_fastq.write(names[i], comments[i], seqs[i], quals[i]);
        }
        writer_fastq.close();
        
        btllib::SeqReader reader_fastq(random_filename + ".bz2");
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

        // Test larger randomly generated file
        std::cerr << "Test random file" << std::endl;
        std::vector<std::string> generated_names;
        std::vector<std::string> generated_comments;
        std::vector<std::string> generated_seqs;
        std::vector<std::string> generated_quals;
        random_filename.clear();
        for (int n = 0; n < 64; n++) {
            random_filename += char(gen_random_alphabet());
        }
        btllib::SeqWriter random_seqs(random_filename, btllib::SeqWriter::FASTQ, false);
        for (int s = 0; s < 500; s++) {
            std::string name, comment, seq, qual;

            for (int n = 0; n < 10; n++) {
                name += char(gen_random_alphabet());
            }
            for (int n = 0; n < 20; n++) {
                comment += char(gen_random_alphabet());
            }
            for (int n = 0; n < 200; n++) {
                seq += "ACTG"[(gen_random_actg())];
            }
            for (int n = 0; n < 200; n++) {
                qual += char(gen_random_alphabet());
            }

            random_seqs.write(name, comment, seq, qual);

            generated_names.push_back(name);
            generated_comments.push_back(comment);
            generated_seqs.push_back(seq);
            generated_quals.push_back(qual);
        }
        random_seqs.close();

        btllib::SeqReader random_reader(random_filename);
        for (int n = 0; record = random_reader.read(); n++) {
            if (record.name != generated_names[n]) { std::cerr << record.name << " | " << generated_names[n] << std::endl; }
            assert(record.name == generated_names[n]);
            assert(record.comment == generated_comments[n]);
            assert(record.seq == generated_seqs[n]);
            assert(record.qual == generated_quals[n]);
        }
    }

	return 0;
}
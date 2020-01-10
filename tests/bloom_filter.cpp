#include "../include/btllib/bloom_filter.hpp"

#include <random>
#include <chrono>
#include <cstdio>

int main() {
    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count() + 7777777);
    std::uniform_int_distribution<int> distribution_alphabet(65, 90);
    auto gen_random_alphabet = std::bind(distribution_alphabet, generator);
    std::string random_filename;

    btllib::BloomFilter bf(1024 * 1024, 3);
    bf.insert({ 1, 10, 100 });
    bf.insert({ 100, 200, 300 });

    assert(bf.contains({ 1, 10, 100 }));
    assert(bf.contains({ 100, 200, 300 }));
    assert(!bf.contains({ 1, 20, 100 }));

    random_filename.clear();
    for (int n = 0; n < 64; n++) {
        random_filename += char(gen_random_alphabet());
    }
    bf.write(random_filename);

    btllib::BloomFilter bf2(random_filename);

    assert(bf2.contains({ 1, 10, 100 }));
    assert(bf2.contains({ 100, 200, 300 }));
    assert(!bf2.contains({ 1, 20, 100 }));

    std::remove(random_filename.c_str());

    return 0;
}
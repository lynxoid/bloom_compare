#include <fstream>
#include <memory>
#include <unordered_set>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <chrono>

// libbf
#include "vallentine/BaseBloomFilter.hpp"
#include "io.hpp"

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
// Usage:
//./sbf <input file type ["reads"|"kmers"]> <input fasta> <query fasta> <filter type ["classic" | "onesided" | "twosided" | "sparse"]> <k> [# queries = 100000]
int main(int argc, char* argv[]) {
    cerr << "===============================" << endl;
    cerr << "Testing Vallentine Bloom Filter" << endl;
    cerr << "===============================" << endl;

    string input_fasta = argv[1];
    int K = stoi(argv[2]);
    cerr << "Counting kmers" << endl;
    shared_ptr<unordered_set<kmer_t>> read_kmers = parseAndCount(input_fasta, 20);
    cerr << "Input kmers: " << read_kmers->size() << endl;

    // allocate more size to avoid collisions
    BaseBloomFilter b(K, read_kmers->size() * 15);
    auto start = std::chrono::system_clock::now();
    // insertions
    for (auto km : *read_kmers)
        b.add(km);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    cerr << "insertion: " << elapsed_seconds.count() << "s" << endl;

    // randomly select 1mln kmers, mutate half of them (FP) and query using this
    // set
    shared_ptr<vector<kmer_t>> query_kmers = select_query_set(read_kmers, 1000000);

    // querying
    start = std::chrono::system_clock::now();
    // queryKmers(query_kmers, read_kmers, b, "classic_query_results.txt");
    for (const kmer_t kmer : *query_kmers)
        b.contains(kmer);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    cerr << "query for 1mln kmers (50% TP, 50% FP): " << elapsed_seconds.count() << "s " << endl;
}

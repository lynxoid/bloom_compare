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
#include <memory>

// libbf
#include "BaseBloomFilter.hpp"
#include "FastaReader.h"
// #include "SeqBFUtil.hpp"
#include "JellyfishUtil.h"

#include <string>

using namespace std;


////////////////////////////////////////////////////////////////////////////////
shared_ptr<unordered_set<kmer_t>> parseAndCount(const string & path,
    const int K) {
    shared_ptr<unordered_set<kmer_t>> kmers(new unordered_set<kmer_t>());
    FastaReader fr(path.c_str());
    kseq_t * seq;
    size_t cnt = 0;
    while ( (seq = fr.nextSequence() ) ) {
      // cerr << seq->seq.s << endl;
        // TODO: INTRODUCED A bug somewhere here
        // string r = seq->seq.s;
        if (seq->seq.l < K) continue;
        for (int i = 0; i < seq->seq.l - K + 1; i++) {
            kmer_t kmer_bin = nimble::mer_string_to_binary(seq->seq.s + i, K);
            kmers->insert(kmer_bin);
        }
        cnt++;
    }
    // cerr << "Parsed " << cnt << " reads" << endl;
    cerr << "(" << cnt << " reads) ";
    return kmers;
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void queryKmers(vector<kmer_t> & test_kmers, unordered_set<kmer_t> & true_kmers,
    BaseBloomFilter & sbf, const string & out_fname) {
    vector<bool> states;
    // time this part
    auto start = std::chrono::system_clock::now();
    for (auto qk : test_kmers) {
        bool state = sbf.contains(qk);
        states.push_back(state);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    cerr << "query: " << elapsed_seconds.count() << "s" << endl;
    // end the timing here

    // write states and true answers to file
    ofstream f_out(out_fname);
    f_out << "kmer\tBF_state\ttrue_state" << endl;
    for (int i = 0; i < states.size(); i++)
        f_out << test_kmers[i] << "\t" <<
          states[i] << "\t" <<
          (true_kmers.find(test_kmers[i]) != true_kmers.end() ) << endl;
    f_out.close();
    // cerr << "Results written to " << out_fname << endl;
}

// simple mutation
kmer_t mutate_kmer(const kmer_t kmer) {
    return !kmer;
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<kmer_t>> select_query_set(
    const shared_ptr<unordered_set<kmer_t>> true_kmers,
    int const set_size) {
    // put mixins into a vector and then randomly draw an index
    shared_ptr<vector<kmer_t>> query_kmers(new vector<kmer_t>());

    // temporarily put into a vector
    vector<kmer_t> temp;
    for (auto k : *true_kmers)
        temp.push_back(k);

    // set up a random number gen
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, true_kmers->size() - 1);
    std::uniform_real_distribution<> mutate_distro(0, 1);
    int i = 0;
    for (int i = 0; i < set_size; i++) {
        auto r = dis(gen);
        assert(r < temp.size());
        // mutate?
        kmer_t kmer = temp[r];
        int result = round(mutate_distro(gen));
        if (result) {
            kmer = mutate_kmer(kmer);
        }
        query_kmers->push_back(kmer);
    }
    return query_kmers;
}

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
    shared_ptr<unordered_set<kmer_t>> read_kmers = parseAndCount(input_fasta, 20);
    cerr << "Input kmers: " << read_kmers->size() << " ";

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

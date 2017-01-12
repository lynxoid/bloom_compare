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
#include "BaseBloomFilter.hpp"
#include "SeqBFUtil.hpp"
#include "JellyfishUtil.h"

#include <string>

using namespace std;


////////////////////////////////////////////////////////////////////////////////
unordered_set<kmer_t> parseAndCount(const string & path,
    const int K) {
    unordered_set<kmer_t> kmers;
    FastaReader fr(path.c_str());
    kseq_t * seq;
    size_t cnt = 0;
    while ( (seq = fr.nextSequence() ) ) {
      // cerr << seq->seq.s << endl;
        // TODO: INTRODUCED A bug somewhere here
        // string r = seq->seq.s;
        if (seq->seq.l < K) continue;
        for (int i = 0; i < seq->seq.l - K + 1; i++) {
            kmer_t kmer_bin = mer_string_to_binary(seq->seq.s + i, K);
            kmers.insert(kmer_bin);
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

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
vector<kmer_t> mix_kmers(unordered_set<kmer_t> & true_kmers, unordered_set<kmer_t> & mixin_kmers,
    int const set_size) {
    // put mixins into a vector and then randomly draw an index
    vector<kmer_t> v;
    for (auto km : mixin_kmers) v.push_back(km);

    vector<kmer_t> query_kmers;

    // set up a random number gen
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, v.size() - 1);
    int i = query_kmers.size();
    while (i < set_size) {
        auto r = dis(gen);
        assert(r < v.size());
        query_kmers.push_back(v[r]);
        i++;
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
    // what's this?
    string mixin_kmers_fasta = argv[2];
    int K = stoi(argv[3]);
    unsigned long query_set_size = 1000000;
    if (argc > 5) {
        query_set_size = stoi(argv[4]);
    }

    unordered_set<kmer_t> read_kmers;
    unordered_set<kmer_t> edge_kmers;
    unordered_set<kmer_t> mixin_kmers;

    ////////////////////////////////////////////////////////
    // parse input reads -- will build a kmer set from that
    ////////////////////////////////////////////////////////
    vector<string> reads = parseAndCount(input_fasta, 20);
    read_kmers = getKmers(reads, K);
    cerr << "Input kmers: " << read_kmers.size() << " ";

    //////////////////////////////////////////////////////
    // test the classic bloom filter (base line)
    //////////////////////////////////////////////////////
    cerr << "Classic: ";
    // allocate more size to avoid collisions
    BaseBloomFilter b(K, read_kmers.size() * 15);
    auto start = std::chrono::system_clock::now();
    // insertions
    for (auto km : read_kmers)
        b.add(km);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    cerr << "build: " << elapsed_seconds.count() << "s ";
    // querying
    queryKmers(query_kmers, read_kmers, b, "classic_query_results.txt");
}

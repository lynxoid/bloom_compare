#ifndef IO_LIB
#define IO_LIB

#include <memory>
#include <fstream>
#include <unordered_set>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <cassert>

#include "FastaReader.h"
// #include "SeqBFUtil.hpp"
#include "JellyfishUtil.h"
#include "definitions.hpp"

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
        if (seq->seq.l < K) continue;
        for (int i = 0; i < seq->seq.l - K + 1; i++) {
            kmer_t kmer_bin = nimble::mer_string_to_binary(seq->seq.s + i, K);
            kmers->insert(kmer_bin);
        }
        cnt++;
    };
    cerr << "(" << cnt << " reads) ";
    return kmers;
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

#endif

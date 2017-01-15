/*
 **************************************************************************
 *                                                                        *
 *                           Open Bloom Filter                            *
 *                                                                        *
 * Description: Basic Bloom Filter Usage                                  *
 * Author: Arash Partow - 2000                                            *
 * URL: http://www.partow.net                                             *
 * URL: http://www.partow.net/programming/hashfunctions/index.html        *
 *                                                                        *
 * Copyright notice:                                                      *
 * Free use of the Bloom Filter Library is permitted under the guidelines *
 * and in accordance with the most current version of the Common Public   *
 * License.                                                               *
 * http://www.opensource.org/licenses/cpl1.0.php                          *
 *                                                                        *
 **************************************************************************
*/



/*
   Description: This example demonstrates basic usage of the Bloom filter.
                Initially some values are inserted then they are subsequently
                queried, noting any false positives or errors.
*/


#include <iostream>
#include <string>
#include <chrono>
#include <vector>
#include <unordered_set>
#include <memory>

#include "partow/bloom_filter.hpp"
#include "definitions.hpp"
#include "io.hpp"

int main(int argc, char * argv[]) {
    cerr << "===============================" << endl;
    cerr << "Testing Partow Bloom Filter" << endl;
    cerr << "===============================" << endl;

    string input_fasta = argv[1];
    int K = stoi(argv[2]);
    cerr << "Counting kmers" << endl;
    shared_ptr<unordered_set<kmer_t>> read_kmers = parseAndCount(input_fasta, 20);
    cerr << "Input kmers: " << read_kmers->size() << endl;

    bloom_parameters parameters;
    // How many elements roughly do we expect to insert?
    parameters.projected_element_count = 10000000;
    // Maximum tolerable false positive probability? (0,1)
    parameters.false_positive_probability = 0.01; // 1 in 100
    // Simple randomizer (optional)
    parameters.random_seed = 0xA5A5A5A5;
    if (!parameters) {
        std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
        return 1;
    }
    parameters.compute_optimal_parameters();
    //Instantiate Bloom Filter
    bloom_filter filter(parameters);

    // Insert into Bloom Filter
    {
        auto start = std::chrono::system_clock::now();
        // Insert some strings
        for (auto kmer : *read_kmers) {
            filter.insert(kmer);
        }
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        cerr << "insertion: " << elapsed_seconds.count() << "s" << endl;
    }

    shared_ptr<vector<kmer_t>> query_kmers = select_query_set(read_kmers, 1000000);
    // Query Bloom Filter
    {
        auto start = std::chrono::system_clock::now();
        // Query the existence of strings
        for (const kmer_t kmer : *query_kmers) {
            filter.contains(kmer);
        }
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        cerr << "query for 1mln kmers (50% TP, 50% FP): " << elapsed_seconds.count() << "s " << endl;
    }
}

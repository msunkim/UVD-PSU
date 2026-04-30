// test/our_verify_sim.cc
//
// Micro-benchmark for our BF-based verification primitive: for each
// recovered union element, query the appropriate Bloom filter (k=64 hash
// inspections per query).  This is the same operation that runs at the
// end of our Decode+Verify phase, isolated and parallelized across 24
// threads to match the wall-clock setting used by prism_verify_sim.cc.
//
// Output is a CSV with the same column layout (n, our_verify_ms_par,
// our_verify_ms_serial) so the two simulations can be plotted/tabled
// side by side.
//
// We populate a single Bloom filter with n random 64-bit values and then
// query it with the same n values (so every query is a true positive,
// matching the intended verification path; the cost of a query is
// independent of membership).

#include "bloom_filter.h"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <omp.h>

int main(int argc, char** argv) {
    std::vector<size_t> ns;
    if (argc > 1) {
        std::string s = argv[1];
        size_t pos = 0;
        while (pos < s.size()) {
            size_t comma = s.find(',', pos);
            ns.push_back(std::stoull(s.substr(pos, comma - pos)));
            if (comma == std::string::npos) break;
            pos = comma + 1;
        }
    } else {
        for (int k = 10; k <= 20; ++k) ns.push_back(size_t{1} << k);
    }

    int trials = 3;
    if (argc > 2) trials = std::atoi(argv[2]);
    const char* csv_path = (argc > 3) ? argv[3] : "our_verify_sim.csv";

    constexpr size_t kHashes = 64;

    std::ofstream csv(csv_path);
    csv << "n,our_verify_ms_par,our_verify_ms_serial\n";

    std::cerr << "Our BF-verification micro-benchmark (k=" << kHashes
              << ", " << omp_get_max_threads() << " threads)\n";
    std::cerr << "n,our_verify_ms_par,our_verify_ms_serial\n";

    for (size_t n : ns) {
        // Build a populated Bloom filter sized for n elements at k=64.
        vdpsu::BloomFilter bf(n, kHashes);

        std::mt19937_64 rng(0xBADCAFEEULL ^ (uint64_t)n);
        std::vector<int64_t> vals(n);
        for (size_t i = 0; i < n; ++i) {
            vals[i] = (int64_t)rng();
            bf.Insert(vals[i]);
        }

        auto run_pass = [&](bool parallel) {
            std::atomic<size_t> hit_counter{0};
            auto t0 = std::chrono::steady_clock::now();
            if (parallel) {
                size_t local_hits = 0;
                #pragma omp parallel reduction(+:local_hits)
                {
                    #pragma omp for schedule(static)
                    for (size_t i = 0; i < n; ++i) {
                        if (bf.Query(vals[i])) ++local_hits;
                    }
                }
                hit_counter.store(local_hits);
            } else {
                size_t hits = 0;
                for (size_t i = 0; i < n; ++i) {
                    if (bf.Query(vals[i])) ++hits;
                }
                hit_counter.store(hits);
            }
            auto t1 = std::chrono::steady_clock::now();
            // hit_counter.load() == n by construction; just consume to
            // prevent the compiler from eliding the calls.
            volatile size_t sink = hit_counter.load();
            (void)sink;
            return std::chrono::duration<double, std::milli>(t1 - t0).count();
        };

        // Warm up
        run_pass(true);

        double sum_par = 0, sum_serial = 0;
        for (int t = 0; t < trials; ++t) sum_par += run_pass(true);
        for (int t = 0; t < trials; ++t) sum_serial += run_pass(false);
        double mean_par = sum_par / trials;
        double mean_serial = sum_serial / trials;

        csv << n << "," << mean_par << "," << mean_serial << "\n";
        csv.flush();
        std::cerr << n << "," << mean_par << "," << mean_serial << "\n";
    }
    return 0;
}

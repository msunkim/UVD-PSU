/**
 * @file bench/main.cc
 * @brief VD-PSU Benchmark: Measures performance of the protocol phases
 *
 * This benchmark faithfully implements the VD-PSU protocol from Section 4
 * and measures timing for each phase: Setup, Encode, Delegate, Compute, Decode.
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <random>
#include <algorithm>
#include <fstream>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "params.h"
#include "prf.h"
#include "mht.h"
#include "phf.h"
#include "timer.h"
#include "fhe_utils.h"
#include "client.h"
#include "cloud.h"

using namespace vdpsu;

// Configuration for benchmark
struct BenchConfig {
    std::vector<size_t> set_sizes;  // n values to test
    size_t num_trials;              // number of repetitions per n
    double intersection_ratio;       // |X_A cap X_B| / n
    bool verbose;
    std::string output_file;
};

// Generate test sets X_A and X_B with specified intersection ratio.
// Elements are 32-bit integers (fit well within P=2^50 plaintext modulus).
// Uses a deterministic contiguous allocation scheme so memory is O(n)
// rather than O(n) in a std::set of ~50M entries (critical at n=2^24).
void GenerateTestSets(size_t n, double intersection_ratio,
                      std::vector<int64_t>& X_A, std::vector<int64_t>& X_B,
                      std::mt19937_64& rng) {
    X_A.clear();
    X_B.clear();

    size_t intersection_size = static_cast<size_t>(n * intersection_ratio);
    size_t unique_per_set = n - intersection_size;

    // Assign disjoint integer ranges that fit in [1, ELEMENT_MAX).
    //   intersection : [1,           1 + intersection_size)
    //   unique_A     : [base_A,      base_A + unique_per_set)
    //   unique_B     : [base_B,      base_B + unique_per_set)
    // Partitioning ELEMENT_MAX into thirds guarantees disjointness for n <= ELEMENT_MAX/3.
    X_A.reserve(n);
    X_B.reserve(n);

    const int64_t third = ELEMENT_MAX / 3;
    const int64_t base_inter = 1;
    const int64_t base_A     = third;
    const int64_t base_B     = 2 * third;
    if (static_cast<int64_t>(std::max(intersection_size, unique_per_set)) >= third) {
        throw std::runtime_error("GenerateTestSets: n too large for ELEMENT_MAX range");
    }

    for (size_t i = 0; i < intersection_size; ++i) {
        int64_t v = base_inter + static_cast<int64_t>(i);
        X_A.push_back(v);
        X_B.push_back(v);
    }
    for (size_t i = 0; i < unique_per_set; ++i) {
        X_A.push_back(base_A + static_cast<int64_t>(i));
        X_B.push_back(base_B + static_cast<int64_t>(i));
    }

    // Shuffle to avoid any ordering bias in PHF / BF insertion
    std::shuffle(X_A.begin(), X_A.end(), rng);
    std::shuffle(X_B.begin(), X_B.end(), rng);
}

// Run a single benchmark trial
BenchmarkMetrics RunTrial(const std::vector<int64_t>& X_A,
                          const std::vector<int64_t>& X_B,
                          FHEUtils& fhe,
                          bool verbose) {
    BenchmarkMetrics metrics;
    Timer timer;
    size_t n = X_A.size();
    ProtocolConfig cfg(n);

    // ========================================
    // Setup: Build shared PHF
    // ========================================
    timer.Start();

    // Create temporary clients to compute witnesses
    PHF shared_phf;
    ClientA clientA(cfg, fhe, shared_phf);
    ClientB clientB(cfg, fhe, shared_phf);

    NTL::ZZ wit_A = clientA.ComputeWitness();
    NTL::ZZ wit_B = clientB.ComputeWitness();

    // Build PHF from union of all elements + witnesses
    std::set<int64_t> all_elems_set;
    for (auto x : X_A) all_elems_set.insert(x);
    for (auto x : X_B) all_elems_set.insert(x);
    all_elems_set.insert(NTL::conv<long>(wit_A));
    all_elems_set.insert(NTL::conv<long>(wit_B));

    std::vector<NTL::ZZ> all_elems;
    for (auto x : all_elems_set) {
        all_elems.push_back(NTL::to_ZZ(static_cast<long>(x)));
    }
    shared_phf.Build(all_elems, cfg.m);

    // Setup time not included in phase metrics (it's a one-time cost)
    double setup_ms = timer.StopMs();
    if (verbose) {
        std::cout << "    Setup (PHF): " << setup_ms << " ms" << std::endl;
    }

    // Create cloud
    Cloud cloud(cfg, fhe);

    // ========================================
    // Phase 1: Encode (Algorithm 1 in paper) — per-client timing
    // In deployment A and B run on separate machines in parallel, so we
    // report each party's local wall-clock time individually.
    // ========================================
    timer.Start();
    auto bar_e_A = clientA.Encode(X_A);
    metrics.encode_A.time_ms = timer.StopMs();

    timer.Start();
    auto bar_e_B = clientB.Encode(X_B);
    metrics.encode_B.time_ms = timer.StopMs();

    // Communication: each client uploads an m-length vector to cloud.
    metrics.encode_A.bytes_sent = cfg.m * sizeof(int64_t);
    metrics.encode_B.bytes_sent = cfg.m * sizeof(int64_t);

    // ========================================
    // Phase 2: Delegate (Algorithm 2 in paper) — per-client timing
    // Note: the generated key pair is passed from A to B, which is a
    // protocol-level message, not a computation; we measure A's KeyGen and
    // Encrypt inside Delegate_A, and B's pure Encrypt inside Delegate_B.
    // ========================================
    timer.Start();
    auto token_A_cloud = clientA.Delegate();
    metrics.delegate_A.time_ms = timer.StopMs();

    timer.Start();
    auto token_B_cloud = clientB.Delegate(clientA.GetKeyPair());
    auto token_B_to_A = clientB.GetTokenForA();
    metrics.delegate_B.time_ms = timer.StopMs();

    size_t ct_size = fhe.CiphertextBytes();
    size_t cts_per_vector = cfg.num_cts;

    // A sends: 2 CT vectors to cloud.
    // (FHE public + relinearization keys are setup payloads, excluded from the per-execution comm.)
    metrics.delegate_A.bytes_sent = 2 * cts_per_vector * ct_size;

    // B sends to A: xmask_bf_B + wit_B (8 bytes), and to cloud: 3 CT vectors.
    const double bits_per_elem = static_cast<double>(BF_NUM_HASHES) / 0.6931471805599453;
    const size_t bf_bits = static_cast<size_t>(bits_per_elem * (n + 1)) + 1;
    const size_t bf_bytes = ((bf_bits + 63) / 64) * 8;
    metrics.delegate_B.bytes_sent = bf_bytes + 8 + 3 * cts_per_vector * ct_size;

    // ========================================
    // Phase 3: Compute (cloud computation)
    // ========================================
    timer.Start();
    auto cloud_response = cloud.Compute(token_A_cloud, token_B_cloud,
                                         bar_e_A, bar_e_B);
    metrics.compute.time_ms = timer.StopMs();
    metrics.compute.bytes_sent = 4 * cts_per_vector * ct_size;

    // ========================================
    // Phase 4: Decode + Verify (Algorithm 3) — client A only
    // ========================================
    timer.Start();
    auto union_result = clientA.Decode(cloud_response, token_B_to_A);
    metrics.decode_verify.time_ms = timer.StopMs();

    if (union_result.empty()) {
        std::cerr << "ERROR: Decode returned empty (verification failure)" << std::endl;
    }

    // ========================================
    // Phase 5: Update (Algorithm 5) — measure mean per-update time on A.
    // We sample up to 100 elements from X_A and time a delete+insert cycle each;
    // the reported update_us is the mean of one (delete or insert) operation.
    // ========================================
    {
        const size_t MAX_SAMPLES = 100;
        size_t samples = std::min<size_t>(MAX_SAMPLES, X_A.size());
        UpdateRequest req;
        std::vector<int64_t> bar_e_A_copy = bar_e_A;  // mutable copy for cloud-side updates
        timer.Start();
        for (size_t k = 0; k < samples; ++k) {
            int64_t elt = X_A[k];
            // Delete then re-insert.
            if (clientA.Update('-', elt, req)) {
                cloud.ProcessUpdate(req, bar_e_A_copy);
            }
            if (clientA.Update('+', elt, req)) {
                cloud.ProcessUpdate(req, bar_e_A_copy);
            }
        }
        double total_ms = timer.StopMs();
        // Two operations per sample (one delete + one insert).
        metrics.update_us = (total_ms * 1000.0) / static_cast<double>(2 * samples);
    }

    // Parallel-deployment wall-clock total:
    //   Encode, Delegate phases are each max(A, B).
    //   Cloud Compute is single-party.
    //   Decode runs only on A.
    auto max_enc = std::max(metrics.encode_A.time_ms, metrics.encode_B.time_ms);
    auto max_del = std::max(metrics.delegate_A.time_ms, metrics.delegate_B.time_ms);
    metrics.total.time_ms = max_enc + max_del + metrics.compute.time_ms +
                            metrics.decode_verify.time_ms;
    metrics.total.bytes_sent = metrics.encode_A.bytes_sent + metrics.encode_B.bytes_sent +
                               metrics.delegate_A.bytes_sent + metrics.delegate_B.bytes_sent +
                               metrics.compute.bytes_sent;

    if (verbose) {
        std::cout << "    Encode_A:      " << std::fixed << std::setprecision(2)
                  << metrics.encode_A.time_ms << " ms" << std::endl;
        std::cout << "    Encode_B:      " << metrics.encode_B.time_ms << " ms" << std::endl;
        std::cout << "    Delegate_A:    " << metrics.delegate_A.time_ms << " ms" << std::endl;
        std::cout << "    Delegate_B:    " << metrics.delegate_B.time_ms << " ms" << std::endl;
        std::cout << "    Compute:       " << metrics.compute.time_ms << " ms" << std::endl;
        std::cout << "    Decode+Verify: " << metrics.decode_verify.time_ms << " ms" << std::endl;
        std::cout << "    Update (mean): " << metrics.update_us << " us/op" << std::endl;
        std::cout << "    Total (wall):  " << metrics.total.time_ms << " ms" << std::endl;
        std::cout << "    Union size:    " << union_result.size() << std::endl;
    }

    return metrics;
}

// Run benchmark for a given set size
void RunBenchmark(size_t n, size_t num_trials, double intersection_ratio,
                  FHEUtils& fhe, bool verbose, std::ostream& out) {
    std::cout << "\n=== Benchmark: n = " << n << " ===" << std::endl;

    std::mt19937_64 rng(42 + n);  // deterministic seed per n for reproducibility

    // Aggregate metrics
    BenchmarkMetrics avg;
    std::vector<BenchmarkMetrics> all_metrics;

    for (size_t trial = 0; trial < num_trials; ++trial) {
        if (verbose) {
            std::cout << "  Trial " << (trial + 1) << "/" << num_trials << ":" << std::endl;
        }

        // Generate test data
        std::vector<int64_t> X_A, X_B;
        GenerateTestSets(n, intersection_ratio, X_A, X_B, rng);

        // Run trial
        auto m = RunTrial(X_A, X_B, fhe, verbose);
        all_metrics.push_back(m);

        avg.encode_A.time_ms     += m.encode_A.time_ms;
        avg.encode_B.time_ms     += m.encode_B.time_ms;
        avg.delegate_A.time_ms   += m.delegate_A.time_ms;
        avg.delegate_B.time_ms   += m.delegate_B.time_ms;
        avg.compute.time_ms      += m.compute.time_ms;
        avg.decode_verify.time_ms+= m.decode_verify.time_ms;
        avg.update_us            += m.update_us;
        avg.total.time_ms        += m.total.time_ms;
    }

    // Compute averages
    avg.encode_A.time_ms     /= num_trials;
    avg.encode_B.time_ms     /= num_trials;
    avg.delegate_A.time_ms   /= num_trials;
    avg.delegate_B.time_ms   /= num_trials;
    avg.compute.time_ms      /= num_trials;
    avg.decode_verify.time_ms/= num_trials;
    avg.update_us            /= num_trials;
    avg.total.time_ms        /= num_trials;

    // Use first trial's communication metrics (they're deterministic)
    avg.encode_A.bytes_sent   = all_metrics[0].encode_A.bytes_sent;
    avg.encode_B.bytes_sent   = all_metrics[0].encode_B.bytes_sent;
    avg.delegate_A.bytes_sent = all_metrics[0].delegate_A.bytes_sent;
    avg.delegate_B.bytes_sent = all_metrics[0].delegate_B.bytes_sent;
    avg.compute.bytes_sent    = all_metrics[0].compute.bytes_sent;
    avg.total.bytes_sent      = all_metrics[0].total.bytes_sent;

    // Print summary
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Average over " << num_trials << " trials:" << std::endl;
    std::cout << "  Encode_A:      " << avg.encode_A.time_ms << " ms" << std::endl;
    std::cout << "  Encode_B:      " << avg.encode_B.time_ms << " ms" << std::endl;
    std::cout << "  Delegate_A:    " << avg.delegate_A.time_ms << " ms" << std::endl;
    std::cout << "  Delegate_B:    " << avg.delegate_B.time_ms << " ms" << std::endl;
    std::cout << "  Compute:       " << avg.compute.time_ms << " ms" << std::endl;
    std::cout << "  Decode+Verify: " << avg.decode_verify.time_ms << " ms" << std::endl;
    std::cout << "  Update (mean): " << avg.update_us << " us/op" << std::endl;
    std::cout << "  Total (wall):  " << avg.total.time_ms << " ms" << std::endl;

    // Communication (in KB/MB) — total across both parties and cloud.
    double total_kb = avg.total.bytes_sent / 1024.0;
    std::cout << "  Total comm:    " << total_kb << " KB";
    if (total_kb > 1024) {
        std::cout << " (" << (total_kb / 1024.0) << " MB)";
    }
    std::cout << std::endl;

    // Output to file/stream (CSV)
    out << n << ","
        << avg.encode_A.time_ms << ","
        << avg.encode_B.time_ms << ","
        << avg.delegate_A.time_ms << ","
        << avg.delegate_B.time_ms << ","
        << avg.compute.time_ms << ","
        << avg.decode_verify.time_ms << ","
        << avg.update_us << ","
        << avg.total.time_ms << ","
        << total_kb << std::endl;
}

void PrintUsage(const char* prog) {
    std::cout << "Usage: " << prog << " [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -n <sizes>      Set sizes to benchmark (comma-separated)" << std::endl;
    std::cout << "                  Default: 1000,2000,4000" << std::endl;
    std::cout << "  -t <trials>     Number of trials per set size (default: 3)" << std::endl;
    std::cout << "  -i <ratio>      Intersection ratio (default: 0.5)" << std::endl;
    std::cout << "  -o <file>       Output CSV file (default: stdout)" << std::endl;
    std::cout << "  -v              Verbose output" << std::endl;
    std::cout << "  -h              Show this help" << std::endl;
}

std::vector<size_t> ParseSizes(const std::string& s) {
    std::vector<size_t> result;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        result.push_back(std::stoull(item));
    }
    return result;
}

int main(int argc, char* argv[]) {
    BenchConfig config;
    config.set_sizes = {1000, 2000, 4000};
    config.num_trials = 3;
    config.intersection_ratio = 0.5;
    config.verbose = false;
    config.output_file = "";

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            PrintUsage(argv[0]);
            return 0;
        } else if (arg == "-n" && i + 1 < argc) {
            config.set_sizes = ParseSizes(argv[++i]);
        } else if (arg == "-t" && i + 1 < argc) {
            config.num_trials = std::stoull(argv[++i]);
        } else if (arg == "-i" && i + 1 < argc) {
            config.intersection_ratio = std::stod(argv[++i]);
        } else if (arg == "-o" && i + 1 < argc) {
            config.output_file = argv[++i];
        } else if (arg == "-v") {
            config.verbose = true;
        }
    }

    std::cout << "==========================================" << std::endl;
    std::cout << "     VD-PSU Protocol Benchmark" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;

    std::cout << "Configuration:" << std::endl;
    std::cout << "  Set sizes: ";
    for (size_t i = 0; i < config.set_sizes.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << config.set_sizes[i];
    }
    std::cout << std::endl;
    std::cout << "  Trials per size: " << config.num_trials << std::endl;
    std::cout << "  Intersection ratio: " << config.intersection_ratio << std::endl;
    std::cout << std::endl;

    std::cout << "BFV Parameters:" << std::endl;
    std::cout << "  Ring dimension N: " << RING_DIM << std::endl;
    std::cout << "  Slots per CT: " << SLOTS_PER_CT << std::endl;
    std::cout << "  Mult depth: " << MULT_DEPTH << std::endl;
    std::cout << "  Plaintext modulus: " << PLAINTEXT_MODULUS
              << " (" << std::fixed << std::setprecision(1)
              << std::log2((double)PLAINTEXT_MODULUS) << " bits)" << std::endl;
#ifdef _OPENMP
    std::cout << "  OpenMP threads: " << omp_get_max_threads() << std::endl;
#endif
    std::cout << std::endl;

    // Initialize FHE context (one-time)
    std::cout << "Initializing FHE context..." << std::endl;
    Timer timer;
    timer.Start();
    FHEUtils fhe;
    fhe.Init();
    double init_ms = timer.StopMs();
    std::cout << "FHE initialization: " << init_ms << " ms" << std::endl;

    // Setup output
    std::ofstream ofs;
    std::ostream* out = &std::cout;
    if (!config.output_file.empty()) {
        ofs.open(config.output_file);
        if (ofs.is_open()) {
            out = &ofs;
            std::cout << "Writing results to: " << config.output_file << std::endl;
        }
    }

    // Write CSV header to output
    *out << "n,encodeA_ms,encodeB_ms,delegateA_ms,delegateB_ms,compute_ms,decode_ms,update_us,total_ms,comm_kb" << std::endl;

    // Run benchmarks
    for (size_t n : config.set_sizes) {
        RunBenchmark(n, config.num_trials, config.intersection_ratio,
                     fhe, config.verbose, *out);
    }

    std::cout << "\n==========================================" << std::endl;
    std::cout << "Benchmark complete." << std::endl;

    return 0;
}

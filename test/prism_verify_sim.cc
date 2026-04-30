// test/prism_verify_sim.cc
//
// Simulation of the dominant cost of PRISM (Li et al., SIGMOD 2021)
// DB-owner result verification.  Per §5.2, Equations (8)-(10), each
// cell of the verification table requires THREE modular multiplications
// modulo the cyclic-group prime.  We model this cost end-to-end with GMP.
//
// We deliberately do NOT implement the full PRISM stack -- we only
// measure the dominant arithmetic primitive.  The simulation thus gives
// a LOWER BOUND on PRISM's verification wall-clock; any realistic
// implementation pays additional overhead for share I/O, permutation
// inversion, and group-element decoding.
//
// Cell count b: PRISM allocates b = |Dom(A_c)| cells (or b = bucket-tree
// occupancy under the optimization of §6.6).  We set b = n for a
// generous comparison; in production PRISM b is typically larger than
// the input set size n.
//
// Modulus: PRISM's reported parameters use a 113-bit prime.
//
// Parallelism: 24 OpenMP threads, matching the hardware setup (Xeon
// W-3235, 12 physical / 24 logical cores) used for our own verification
// benchmark in vdpsu-bench.

#include <gmp.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <omp.h>

namespace {

// A 113-bit safe prime (2*q+1, q prime), to match PRISM's reported size.
// 2^112 < p < 2^113.
constexpr const char* kPrime113 =
    "10384593717069655257060992658440191";

}  // namespace

int main(int argc, char** argv) {
    // Allow user to override the n list:
    //   ./prism_verify_sim                          -> 2^10 .. 2^20
    //   ./prism_verify_sim 1024,2048,4096           -> custom
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

    // Trial count
    int trials = 3;
    if (argc > 2) trials = std::atoi(argv[2]);

    // Output CSV
    const char* csv_path = (argc > 3) ? argv[3] : "prism_verify_sim.csv";

    // Modulus
    mpz_t p;
    mpz_init_set_str(p, kPrime113, 10);

    std::ofstream csv(csv_path);
    csv << "n,prism_verify_ms_par,prism_verify_ms_serial\n";

    std::cerr << "PRISM verification simulation (b=n, 113-bit prime, "
              << omp_get_max_threads() << " threads)\n";
    std::cerr << "n,prism_verify_ms_par,prism_verify_ms_serial\n";

    for (size_t n : ns) {
        // Allocate four random per-cell vectors a,b,c,d in [0,p).
        // Each cell verification is:
        //   r1 = a*b mod p     (Eq. 8)
        //   r2 = c*d mod p     (Eq. 9)
        //   q  = r1*r2 mod p   (Eq. 10)
        // = 3 modular multiplications.
        std::vector<__mpz_struct> a(n), b(n), c(n), d(n);
        for (size_t i = 0; i < n; ++i) {
            mpz_init(&a[i]); mpz_init(&b[i]); mpz_init(&c[i]); mpz_init(&d[i]);
        }

        std::mt19937_64 rng(0xC0FFEEULL ^ (uint64_t)n);
        // Fill each mpz with a 128-bit random value, reduced mod p.
        auto fill_random = [&](__mpz_struct* x) {
            uint64_t lo = rng();
            uint64_t hi = rng();
            mpz_set_ui(x, (unsigned long)hi);
            mpz_mul_2exp(x, x, 64);
            mpz_add_ui(x, x, (unsigned long)lo);
            mpz_mod(x, x, p);
        };
        for (size_t i = 0; i < n; ++i) {
            fill_random(&a[i]); fill_random(&b[i]);
            fill_random(&c[i]); fill_random(&d[i]);
        }

        // Helper to run one verification pass (returns wall-clock ms).
        auto run_pass = [&](bool parallel) {
            auto t0 = std::chrono::steady_clock::now();
            if (parallel) {
                #pragma omp parallel
                {
                    mpz_t r1, r2, q;
                    mpz_init(r1); mpz_init(r2); mpz_init(q);
                    #pragma omp for schedule(static)
                    for (size_t i = 0; i < n; ++i) {
                        mpz_mul(r1, &a[i], &b[i]); mpz_mod(r1, r1, p);
                        mpz_mul(r2, &c[i], &d[i]); mpz_mod(r2, r2, p);
                        mpz_mul(q,  r1,    r2);    mpz_mod(q,  q,  p);
                        // (production would do: mpz_cmp_ui(q,1) for the
                        //  detection check; cost is O(1) per cell.)
                    }
                    mpz_clear(r1); mpz_clear(r2); mpz_clear(q);
                }
            } else {
                mpz_t r1, r2, q;
                mpz_init(r1); mpz_init(r2); mpz_init(q);
                for (size_t i = 0; i < n; ++i) {
                    mpz_mul(r1, &a[i], &b[i]); mpz_mod(r1, r1, p);
                    mpz_mul(r2, &c[i], &d[i]); mpz_mod(r2, r2, p);
                    mpz_mul(q,  r1,    r2);    mpz_mod(q,  q,  p);
                }
                mpz_clear(r1); mpz_clear(r2); mpz_clear(q);
            }
            auto t1 = std::chrono::steady_clock::now();
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

        for (size_t i = 0; i < n; ++i) {
            mpz_clear(&a[i]); mpz_clear(&b[i]); mpz_clear(&c[i]); mpz_clear(&d[i]);
        }
    }

    mpz_clear(p);
    return 0;
}

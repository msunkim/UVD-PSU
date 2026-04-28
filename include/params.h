#ifndef VDPSU_PARAMS_H
#define VDPSU_PARAMS_H

#include <cstdint>
#include <cstddef>

namespace vdpsu {

// Security parameter for PRF (bits)
constexpr size_t SECURITY_PARAM = 256;

// BFV parameters
constexpr size_t RING_DIM = 16384;         // N = 16384 (required for P ~ 2^59 at 128-bit security)
constexpr size_t SLOTS_PER_CT = RING_DIM;  // full N slots with batching (16384)
constexpr size_t MULT_DEPTH = 2;           // depth-2: ct-ct mult (depth 1) + pt-ct mult (depth 2)

// Plaintext modulus: depends on the OpenFHE NATIVEINT build.
//
// VDPSU_NATIVE128 (set at build time) switches between two configurations:
//
//   NATIVEINT=64  (default): P = 50-bit prime, 32-bit elements,
//                 ~17-bit statistical slack. Fast path; n up to 2^24.
//
//   NATIVEINT=128: P = 54-bit prime, 22-bit elements,
//                 ~31-bit statistical slack (54 - 22 - 1).
//                 Note: empirical sweep shows that BFV depth-2 ct*ct
//                 multiplication breaks above 54-bit P even with BEHZ and
//                 NATIVEINT=128 (noise budget exhausted). 54 bits is the
//                 practical ceiling, not the 60-bit RNS MAX_SIZE.
//                 22-bit elements give 4x headroom over n=2^20 = 1,048,576.
//
// Both primes satisfy P ≡ 1 mod 2N = 1 mod 32768 for N=16384 SIMD batching.
// Both fit in signed int64_t, so client arithmetic can stay in int64_t / __int128.
#ifdef VDPSU_NATIVE128
constexpr uint64_t PLAINTEXT_MODULUS = 1125899908022273ULL;     // 2^51 prime (stable in NATIVEINT=128)
constexpr int64_t  ELEMENT_MAX       = (1LL << 22);             // 22-bit elements
#else
constexpr uint64_t PLAINTEXT_MODULUS = 1125899908022273ULL;     // 2^50 prime
constexpr int64_t  ELEMENT_MAX       = (1LL << 32);             // 32-bit elements
#endif

// Bloom filter parameters
constexpr size_t BF_NUM_HASHES = 64;       // k = 64 hash functions

// MHT size factor: m = MHT_FACTOR * n
constexpr size_t MHT_FACTOR = 2;           // m = 2n

// Protocol configuration
struct ProtocolConfig {
    size_t n;           // set size
    size_t m;           // MHT size (= MHT_FACTOR * n)
    size_t k;           // number of BF hash functions
    size_t num_cts;     // number of ciphertexts per MHT vector (= ceil(m / SLOTS_PER_CT))

    explicit ProtocolConfig(size_t set_size)
        : n(set_size),
          m(MHT_FACTOR * set_size),
          k(BF_NUM_HASHES),
          num_cts((MHT_FACTOR * set_size + SLOTS_PER_CT - 1) / SLOTS_PER_CT)
    {}
};

// Timing metrics per phase
struct PhaseMetrics {
    double time_ms = 0.0;
    size_t bytes_sent = 0;
    size_t bytes_recv = 0;
};

struct BenchmarkMetrics {
    // Per-client timings for the phases where both clients do work independently.
    PhaseMetrics encode_A;
    PhaseMetrics encode_B;
    PhaseMetrics delegate_A;
    PhaseMetrics delegate_B;
    // Cloud and decode are single-party phases.
    PhaseMetrics compute;
    PhaseMetrics decode_verify;
    // Mean per-update wall-clock time (Algorithm 5), microseconds.
    double update_us = 0.0;
    PhaseMetrics total;
};

} // namespace vdpsu

#endif // VDPSU_PARAMS_H

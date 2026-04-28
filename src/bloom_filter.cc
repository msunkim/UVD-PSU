#include "bloom_filter.h"

#include <cmath>

namespace vdpsu {

BloomFilter::BloomFilter(size_t expected_elements, size_t num_hashes)
    : k_(num_hashes) {
    // Optimal bits-per-element for k hashes: k / ln(2) ≈ 1.4427 * k
    const double bits_per_elem = static_cast<double>(k_) / 0.6931471805599453;
    size_t bits = static_cast<size_t>(bits_per_elem * expected_elements) + 1;
    // Minimum size guard: at least one 64-bit word
    if (bits < 64) bits = 64;
    // Round up to multiple of 64 bits for word-aligned storage
    size_t words = (bits + 63) / 64;
    bits_.assign(words, 0);
    num_bits_ = words * 64;
}

// Lemire's fast range mapping: given h uniform in [0, 2^64), produces
// a value uniform in [0, n) using one 64x64->128 multiply and a shift.
// Much faster than `h % n` (no division), no power-of-two constraint.
static inline uint64_t FastRange(uint64_t h, uint64_t n) {
    return static_cast<uint64_t>((static_cast<__uint128_t>(h) * n) >> 64);
}

void BloomFilter::Insert(int64_t val) {
    if (num_bits_ == 0) return;
    const uint64_t x = static_cast<uint64_t>(val);
    const uint64_t h1 = SplitMix64(x);
    const uint64_t h2 = SplitMix64(x ^ 0xD1B54A32D192ED03ULL) | 1ULL;  // keep h2 odd
    uint64_t h = h1;
    for (size_t i = 0; i < k_; ++i) {
        const uint64_t idx = FastRange(h, num_bits_);
        // Atomic fetch-or so concurrent OpenMP inserts into the same BF are safe.
        // __ATOMIC_RELAXED is sufficient: BF only cares that every bit-set eventually lands.
        __atomic_fetch_or(&bits_[idx >> 6], (1ULL << (idx & 63)), __ATOMIC_RELAXED);
        h += h2;
    }
}

bool BloomFilter::Query(int64_t val) const {
    if (num_bits_ == 0) return false;
    const uint64_t x = static_cast<uint64_t>(val);
    const uint64_t h1 = SplitMix64(x);
    const uint64_t h2 = SplitMix64(x ^ 0xD1B54A32D192ED03ULL) | 1ULL;
    uint64_t h = h1;
    for (size_t i = 0; i < k_; ++i) {
        const uint64_t idx = FastRange(h, num_bits_);
        if ((bits_[idx >> 6] & (1ULL << (idx & 63))) == 0) return false;
        h += h2;
    }
    return true;
}

}  // namespace vdpsu

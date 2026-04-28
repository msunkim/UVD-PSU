#ifndef VDPSU_BLOOM_FILTER_H
#define VDPSU_BLOOM_FILTER_H

#include <cstdint>
#include <cstddef>
#include <vector>

namespace vdpsu {

// Bit-array Bloom filter with double-hashing (Kirsch-Mitzenmacher).
// Hash function: two 64-bit SplitMix64 mixes of the input value, combined
// as h_i(x) = (h1(x) + i*h2(x)) mod m.
//
// For k hash functions, the optimal bit count is m ≈ 1.4427 * k * n where
// n is the expected number of inserted elements. The false positive rate at
// full load is then ≈ 2^-k.
//
// Copyable so it can be passed by value inside TokenB_ToA.
class BloomFilter {
public:
    BloomFilter() = default;
    BloomFilter(size_t expected_elements, size_t num_hashes);

    void Insert(int64_t val);
    bool Query(int64_t val) const;

    size_t NumBits() const { return num_bits_; }
    size_t NumHashes() const { return k_; }
    size_t ByteSize() const { return bits_.size() * sizeof(uint64_t); }

private:
    size_t num_bits_ = 0;
    size_t k_ = 0;
    std::vector<uint64_t> bits_;  // packed bit array

    static inline uint64_t SplitMix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
};

}  // namespace vdpsu

#endif  // VDPSU_BLOOM_FILTER_H

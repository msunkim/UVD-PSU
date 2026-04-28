#ifndef VDPSU_PHF_H
#define VDPSU_PHF_H

#include <cstdint>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <NTL/ZZ.h>

namespace vdpsu {

// Perfect Hash Function (PHF)
// Maps a set of n elements injectively to [0..m-1] where m >= 2n.
// For benchmarking, we use a simple approach: store the mapping explicitly.
// In production, one would use CMPH (CHD/BDZ algorithm).
class PHF {
public:
    PHF() = default;

    // Build PHF from a set of elements, mapping to [0..m-1]
    // Elements are represented as string keys (hex of ZZ)
    void Build(const std::vector<NTL::ZZ>& elements, size_t m);

    // Evaluate: PHF(x) -> index in [0..m-1]
    size_t Eval(const NTL::ZZ& x) const;

    // Fast path for int64_t elements: skips NTL::ZZ construction and string
    // serialization. Requires the element to have been inserted during Build.
    size_t Eval(int64_t x) const;

    // Check if element is in the domain
    bool Contains(const NTL::ZZ& x) const;

    // Get the inverse mapping: index -> element (for decoding)
    const NTL::ZZ& InverseEval(size_t index) const;

    // Get the domain size
    size_t DomainSize() const { return forward_.size(); }

    // Get the range size
    size_t RangeSize() const { return m_; }

private:
    size_t m_ = 0;
    std::unordered_map<std::string, size_t> forward_;     // element -> index (by string)
    std::unordered_map<int64_t, size_t> forward_int_;     // fast path for int64 elements
    std::vector<NTL::ZZ> inverse_;                        // index -> element (sparse)
    std::vector<bool> used_;                              // which indices are occupied

    static std::string ZZToKey(const NTL::ZZ& z);
};

} // namespace vdpsu

#endif // VDPSU_PHF_H

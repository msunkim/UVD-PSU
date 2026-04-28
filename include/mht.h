#ifndef VDPSU_MHT_H
#define VDPSU_MHT_H

#include <cstdint>
#include <cstddef>
#include <vector>
#include <functional>

namespace vdpsu {

// Multiplicity Hash Table (MHT)
// A vector of length m storing multiplicity counts (0 or 1 for sets).
// Indexed by a Perfect Hash Function (PHF).
class MHT {
public:
    // Create an MHT of size m, initialized to all zeros
    explicit MHT(size_t m);

    // Insert element x using hash function: mht[hash(x)] = 1
    void Insert(size_t index);

    // Remove element at index: mht[index] = 0
    void Remove(size_t index);

    // Query: return mht[index]
    int64_t Query(size_t index) const;

    // Get the raw vector (for encryption)
    const std::vector<int64_t>& GetVector() const { return data_; }

    // Get a mutable reference (for building tilde versions)
    std::vector<int64_t>& GetMutableVector() { return data_; }

    // Size
    size_t Size() const { return data_.size(); }

    // Compute component-wise product with another vector: result[j] = data_[j] * other[j]
    std::vector<int64_t> HadamardProduct(const std::vector<int64_t>& other) const;

    // Compute component-wise addition with another MHT
    MHT operator+(const MHT& other) const;

    // Negate: multiply all entries by -1
    MHT Negate() const;

private:
    std::vector<int64_t> data_;
};

} // namespace vdpsu

#endif // VDPSU_MHT_H

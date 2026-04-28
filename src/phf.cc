#include "phf.h"
#include <algorithm>
#include <random>
#include <sstream>
#include <stdexcept>

namespace vdpsu {

std::string PHF::ZZToKey(const NTL::ZZ& z) {
    std::ostringstream oss;
    oss << z;
    return oss.str();
}

void PHF::Build(const std::vector<NTL::ZZ>& elements, size_t m) {
    if (elements.size() > m) {
        throw std::invalid_argument("PHF::Build: more elements than range size");
    }

    m_ = m;
    forward_.clear();
    forward_int_.clear();
    forward_int_.reserve(elements.size());
    inverse_.resize(m);
    used_.assign(m, false);

    // Simple injective mapping: assign each element to a unique random slot
    // For production, use CMPH. For benchmarking, this explicit approach works.
    std::vector<size_t> slots(m);
    std::iota(slots.begin(), slots.end(), 0);

    // Shuffle to get random injective mapping
    std::mt19937_64 rng(42);  // deterministic for reproducibility
    std::shuffle(slots.begin(), slots.end(), rng);

    for (size_t i = 0; i < elements.size(); ++i) {
        size_t slot = slots[i];
        std::string key = ZZToKey(elements[i]);
        forward_[key] = slot;
        inverse_[slot] = elements[i];
        used_[slot] = true;
        // Populate fast int64 path if the element fits in a signed 64-bit integer
        if (NTL::NumBits(elements[i]) <= 63) {
            int64_t iv = NTL::conv<long>(elements[i]);
            forward_int_[iv] = slot;
        }
    }
}

size_t PHF::Eval(const NTL::ZZ& x) const {
    // Fast path: if x fits in int64_t and forward_int_ is populated, use it.
    if (!forward_int_.empty() && NTL::NumBits(x) <= 63) {
        int64_t iv = NTL::conv<long>(x);
        auto it = forward_int_.find(iv);
        if (it != forward_int_.end()) return it->second;
    }
    // Fallback: string-keyed path.
    std::string key = ZZToKey(x);
    auto it = forward_.find(key);
    if (it == forward_.end()) {
        throw std::runtime_error("PHF::Eval: element not in domain");
    }
    return it->second;
}

size_t PHF::Eval(int64_t x) const {
    auto it = forward_int_.find(x);
    if (it == forward_int_.end()) {
        throw std::runtime_error("PHF::Eval(int64_t): element not in domain");
    }
    return it->second;
}

bool PHF::Contains(const NTL::ZZ& x) const {
    std::string key = ZZToKey(x);
    return forward_.find(key) != forward_.end();
}

const NTL::ZZ& PHF::InverseEval(size_t index) const {
    if (index >= m_ || !used_[index]) {
        throw std::out_of_range("PHF::InverseEval: index not occupied");
    }
    return inverse_[index];
}

} // namespace vdpsu

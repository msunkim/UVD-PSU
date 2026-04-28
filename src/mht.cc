#include "mht.h"
#include <stdexcept>

namespace vdpsu {

MHT::MHT(size_t m) : data_(m, 0) {}

void MHT::Insert(size_t index) {
    if (index >= data_.size()) {
        throw std::out_of_range("MHT::Insert: index out of range");
    }
    data_[index] = 1;
}

void MHT::Remove(size_t index) {
    if (index >= data_.size()) {
        throw std::out_of_range("MHT::Remove: index out of range");
    }
    data_[index] = 0;
}

int64_t MHT::Query(size_t index) const {
    if (index >= data_.size()) {
        throw std::out_of_range("MHT::Query: index out of range");
    }
    return data_[index];
}

std::vector<int64_t> MHT::HadamardProduct(const std::vector<int64_t>& other) const {
    if (other.size() != data_.size()) {
        throw std::invalid_argument("MHT::HadamardProduct: size mismatch");
    }
    std::vector<int64_t> result(data_.size());
    for (size_t i = 0; i < data_.size(); ++i) {
        result[i] = data_[i] * other[i];
    }
    return result;
}

MHT MHT::operator+(const MHT& other) const {
    if (other.data_.size() != data_.size()) {
        throw std::invalid_argument("MHT::operator+: size mismatch");
    }
    MHT result(data_.size());
    for (size_t i = 0; i < data_.size(); ++i) {
        result.data_[i] = data_[i] + other.data_[i];
    }
    return result;
}

MHT MHT::Negate() const {
    MHT result(data_.size());
    for (size_t i = 0; i < data_.size(); ++i) {
        result.data_[i] = -data_[i];
    }
    return result;
}

} // namespace vdpsu

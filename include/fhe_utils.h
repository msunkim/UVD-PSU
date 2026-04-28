#ifndef VDPSU_FHE_UTILS_H
#define VDPSU_FHE_UTILS_H

#include "params.h"
#include <vector>
#include <cstdint>

// OpenFHE headers
#include "openfhe.h"

namespace vdpsu {

using CC = lbcrypto::CryptoContext<lbcrypto::DCRTPoly>;
using KeyPair = lbcrypto::KeyPair<lbcrypto::DCRTPoly>;
using CT = lbcrypto::Ciphertext<lbcrypto::DCRTPoly>;
using PT = lbcrypto::Plaintext;

// A packed ciphertext vector: one MHT of length m is stored
// as ceil(m / SLOTS_PER_CT) ciphertexts.
using CTVector = std::vector<CT>;

// BFV FHE utilities for VD-PSU
class FHEUtils {
public:
    FHEUtils() = default;

    // Initialize BFV context with SIMD batching
    // plaintext_modulus: must be prime and = 1 mod 2N for batching
    // If 0, a suitable default is chosen
    void Init(uint64_t plaintext_modulus = 0);

    // Generate key pair (includes relinearization keys for ct-ct multiplication)
    KeyPair KeyGen();

    // --- Packing / Unpacking ---

    // Pack a vector of int64_t (length m) into ceil(m / SLOTS_PER_CT) ciphertexts
    CTVector Encrypt(const KeyPair& kp, const std::vector<int64_t>& vec) const;

    // Decrypt a packed ciphertext vector back to a vector of int64_t
    std::vector<int64_t> Decrypt(const KeyPair& kp, const CTVector& ct_vec, size_t m) const;

    // Encode a plaintext vector (for plaintext-ciphertext operations)
    std::vector<PT> EncodePlaintext(const std::vector<int64_t>& vec) const;

    // --- Homomorphic Operations on Packed Vectors ---

    // Component-wise addition: ct_a[j] + ct_b[j]
    CTVector HAdd(const CTVector& a, const CTVector& b) const;

    // Component-wise ciphertext-ciphertext multiplication: ct_a[j] * ct_b[j] (depth-1)
    CTVector HMul(const CTVector& a, const CTVector& b) const;

    // Component-wise plaintext-ciphertext multiplication: pt[j] * ct[j]
    CTVector PTMul(const std::vector<int64_t>& pt_vec, const CTVector& ct_vec) const;

    // Scalar-ciphertext multiplication: scalar * ct[j] for all j
    CTVector ScalarMul(int64_t scalar, const CTVector& ct_vec) const;

    // --- Utility ---

    // Get the number of ciphertexts needed for a vector of length m
    size_t NumCiphertexts(size_t m) const;

    // Get estimated ciphertext size in bytes
    size_t CiphertextBytes() const;

    // Get the crypto context (for advanced operations)
    CC GetContext() const { return cc_; }

    // Diagnostics
    size_t GetActualRingDim() const;
    size_t GetSlotsPerCT() const;

private:
    CC cc_ = nullptr;
    size_t slots_per_ct_ = SLOTS_PER_CT;
    size_t actual_ring_dim_ = 0;

    // Split a long vector into chunks of size slots_per_ct_
    std::vector<std::vector<int64_t>> SplitVector(
        const std::vector<int64_t>& vec) const;
};

} // namespace vdpsu

#endif // VDPSU_FHE_UTILS_H

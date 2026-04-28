#include "fhe_utils.h"
#include <stdexcept>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

// OpenFHE NATIVEINT=128 (at least in 1.5.0) is NOT thread-safe for concurrent
// EvalAdd / EvalMult / Encrypt / Decrypt on the same crypto context — it
// produces intermittent ciphertext corruption when multiple cts are processed
// in parallel. We serialize the outer FHE loops in that build.
// NATIVEINT=64 has no such issue.
#if defined(VDPSU_NATIVE128)
    #define VDPSU_FHE_OMP_FOR
#else
    #define VDPSU_FHE_OMP_FOR _Pragma("omp parallel for")
#endif

namespace vdpsu {

void FHEUtils::Init(uint64_t plaintext_modulus) {
    // BFV parameters for depth-1 multiplication with SIMD batching
    lbcrypto::CCParams<lbcrypto::CryptoContextBFVRNS> params;

    params.SetRingDim(RING_DIM);
    params.SetMultiplicativeDepth(MULT_DEPTH);
    params.SetSecurityLevel(lbcrypto::HEStd_128_classic);

#ifdef VDPSU_NATIVE128
    // NATIVEINT=128 mode: use the default HPSPOVERQLEVELED multiplication
    // technique, which is more stable at higher num_cts than BEHZ.
    // Constrain plaintext modulus to 51 bits (same as 64-bit build) for stability.
    params.SetScalingModSize(60);
#endif

    // Plaintext modulus: 59-bit prime = 1 mod 2N = 1 mod 32768 for N=16384.
    // Supports 32-bit elements with ~27-bit statistical slack.
    if (plaintext_modulus == 0) {
        plaintext_modulus = PLAINTEXT_MODULUS;
    }
    params.SetPlaintextModulus(plaintext_modulus);

    // Generate crypto context
    cc_ = lbcrypto::GenCryptoContext(params);
    cc_->Enable(lbcrypto::PKE);
    cc_->Enable(lbcrypto::KEYSWITCH);
    cc_->Enable(lbcrypto::LEVELEDSHE);

    // Verify batching is enabled
    auto encodingParams = cc_->GetEncodingParams();
    size_t batchSize = encodingParams->GetBatchSize();
    if (batchSize == 0) {
        throw std::runtime_error("FHEUtils::Init: SIMD batching not enabled. "
                                 "Plaintext modulus must be prime and = 1 mod 2N.");
    }
    slots_per_ct_ = batchSize;

    // Report the actual ring dimension chosen by OpenFHE
    actual_ring_dim_ = cc_->GetCryptoParameters()->GetElementParams()->GetRingDimension();
}

size_t FHEUtils::GetActualRingDim() const {
    return actual_ring_dim_;
}

size_t FHEUtils::GetSlotsPerCT() const {
    return slots_per_ct_;
}

KeyPair FHEUtils::KeyGen() {
    if (!cc_) {
        throw std::runtime_error("FHEUtils::KeyGen: context not initialized");
    }
    auto kp = cc_->KeyGen();
    // Generate relinearization keys for ciphertext-ciphertext multiplication
    cc_->EvalMultKeyGen(kp.secretKey);
    return kp;
}

// --- Packing / Unpacking ---

std::vector<std::vector<int64_t>> FHEUtils::SplitVector(
    const std::vector<int64_t>& vec) const {

    std::vector<std::vector<int64_t>> chunks;
    size_t n = vec.size();
    for (size_t offset = 0; offset < n; offset += slots_per_ct_) {
        size_t chunk_size = std::min(slots_per_ct_, n - offset);
        std::vector<int64_t> chunk(vec.begin() + offset,
                                    vec.begin() + offset + chunk_size);
        // Pad to full slot count if needed
        chunk.resize(slots_per_ct_, 0);
        chunks.push_back(std::move(chunk));
    }
    return chunks;
}

CTVector FHEUtils::Encrypt(const KeyPair& kp,
                           const std::vector<int64_t>& vec) const {
    auto chunks = SplitVector(vec);
    CTVector ct_vec(chunks.size());

    VDPSU_FHE_OMP_FOR
    for (size_t i = 0; i < chunks.size(); ++i) {
        PT pt = cc_->MakePackedPlaintext(chunks[i]);
        ct_vec[i] = cc_->Encrypt(kp.publicKey, pt);
    }
    return ct_vec;
}

std::vector<int64_t> FHEUtils::Decrypt(const KeyPair& kp,
                                        const CTVector& ct_vec,
                                        size_t m) const {
    // Decrypt each ciphertext in parallel into a temporary buffer
    std::vector<std::vector<int64_t>> chunks(ct_vec.size());

    VDPSU_FHE_OMP_FOR
    for (size_t i = 0; i < ct_vec.size(); ++i) {
        PT pt;
        cc_->Decrypt(kp.secretKey, ct_vec[i], &pt);
        pt->SetLength(slots_per_ct_);
        chunks[i] = pt->GetPackedValue();
    }

    // Concatenate (sequential, cheap)
    std::vector<int64_t> result;
    result.reserve(m);
    for (const auto& chunk : chunks) {
        result.insert(result.end(), chunk.begin(), chunk.end());
    }
    result.resize(m);
    return result;
}

std::vector<PT> FHEUtils::EncodePlaintext(const std::vector<int64_t>& vec) const {
    auto chunks = SplitVector(vec);
    std::vector<PT> pt_vec;
    pt_vec.reserve(chunks.size());
    for (const auto& chunk : chunks) {
        pt_vec.push_back(cc_->MakePackedPlaintext(chunk));
    }
    return pt_vec;
}

// --- Homomorphic Operations ---

CTVector FHEUtils::HAdd(const CTVector& a, const CTVector& b) const {
    if (a.size() != b.size()) {
        throw std::invalid_argument("FHEUtils::HAdd: vector size mismatch");
    }
    CTVector result(a.size());
    VDPSU_FHE_OMP_FOR
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = cc_->EvalAdd(a[i], b[i]);
    }
    return result;
}

CTVector FHEUtils::HMul(const CTVector& a, const CTVector& b) const {
    if (a.size() != b.size()) {
        throw std::invalid_argument("FHEUtils::HMul: vector size mismatch");
    }
    CTVector result(a.size());
    VDPSU_FHE_OMP_FOR
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = cc_->EvalMult(a[i], b[i]);
    }
    return result;
}

CTVector FHEUtils::PTMul(const std::vector<int64_t>& pt_vec,
                          const CTVector& ct_vec) const {
    auto pt_encoded = EncodePlaintext(pt_vec);
    if (pt_encoded.size() != ct_vec.size()) {
        throw std::invalid_argument("FHEUtils::PTMul: vector size mismatch");
    }
    CTVector result(ct_vec.size());
    VDPSU_FHE_OMP_FOR
    for (size_t i = 0; i < ct_vec.size(); ++i) {
        result[i] = cc_->EvalMult(ct_vec[i], pt_encoded[i]);
    }
    return result;
}

CTVector FHEUtils::ScalarMul(int64_t scalar, const CTVector& ct_vec) const {
    // Create a plaintext with scalar in all slots
    std::vector<int64_t> scalar_vec(slots_per_ct_, scalar);
    PT pt = cc_->MakePackedPlaintext(scalar_vec);

    CTVector result(ct_vec.size());
    VDPSU_FHE_OMP_FOR
    for (size_t i = 0; i < ct_vec.size(); ++i) {
        result[i] = cc_->EvalMult(ct_vec[i], pt);
    }
    return result;
}

// --- Utility ---

size_t FHEUtils::NumCiphertexts(size_t m) const {
    return (m + slots_per_ct_ - 1) / slots_per_ct_;
}

size_t FHEUtils::CiphertextBytes() const {
    // Approximate: each ciphertext is 2 ring elements of size N * log2(q) bits
    // For depth-1 BFV with N=8192 and ~109-bit modulus: ~2 * 8192 * 14 bytes ≈ 224 KB
    // This is a rough estimate; actual size depends on modulus chain
    return 2 * RING_DIM * 14;  // ~229,376 bytes per ciphertext
}

} // namespace vdpsu

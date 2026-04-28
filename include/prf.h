#ifndef VDPSU_PRF_H
#define VDPSU_PRF_H

#include <cstdint>
#include <vector>
#include <string>
#include <NTL/ZZ.h>

namespace vdpsu {

// PRF key types matching the protocol
enum class KeyType { K1 = 0, K2 = 1, K3 = 2 };

// PRF based on HMAC-SHA256 (Crypto++)
// Evaluates F(K, index) -> 256-bit output in NTL::ZZ
class PRF {
public:
    PRF();
    ~PRF() = default;

    // Generate three PRF keys (K1, K2, K3) from a random seed
    void KeyGen();

    // Generate keys from a specific seed (for reproducibility)
    void KeyGen(const NTL::ZZ& seed);

    // Evaluate PRF: F(key_type, index) -> ZZ
    NTL::ZZ Eval(KeyType type, int32_t index) const;

    // Evaluate PRF and reduce modulo p: F(key_type, index) mod p
    NTL::ZZ EvalMod(KeyType type, int32_t index, const NTL::ZZ& p) const;

    // Fast 64-bit path: returns F(key_type, index) mod p as int64_t.
    // Avoids all NTL::ZZ allocation and uses cached key bytes.
    // Uses the first 8 bytes of the HMAC-SHA256 output; statistical bias
    // from (2^64 mod p) / 2^64 is negligible for p >> 1.
    int64_t EvalMod64(KeyType type, int32_t index, int64_t p) const;

    // Get key (for debugging/testing)
    const NTL::ZZ& GetKey(KeyType type) const;

private:
    NTL::ZZ keys_[3];  // K1, K2, K3
    unsigned char key_bytes_[3][32] = {};  // cached key byte serialization
    bool initialized_ = false;

    // Core HMAC-SHA256 evaluation
    NTL::ZZ HmacEval(const NTL::ZZ& key, int32_t index) const;
};

} // namespace vdpsu

#endif // VDPSU_PRF_H

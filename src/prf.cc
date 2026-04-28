#include "prf.h"

#include <cryptopp/hmac.h>
#include <cryptopp/sha.h>
#include <NTL/ZZ.h>
#include <cstring>
#include <stdexcept>

namespace vdpsu {

static constexpr size_t KEY_BYTES = 32;   // 256-bit key
static constexpr size_t MSG_BYTES = 4;    // 32-bit index

PRF::PRF() : initialized_(false) {}

void PRF::KeyGen() {
    NTL::ZZ seed;
    NTL::RandomLen(seed, 256);
    KeyGen(seed);
}

void PRF::KeyGen(const NTL::ZZ& seed) {
    NTL::ZZ s = seed;
    for (int i = 0; i < 3; ++i) {
        NTL::SetSeed(s);
        NTL::RandomLen(keys_[i], 256);
        s += NTL::to_ZZ(1);
    }
    // Cache the serialized key bytes for the fast int64 path.
    for (int i = 0; i < 3; ++i) {
        std::memset(key_bytes_[i], 0, 32);
        NTL::BytesFromZZ(key_bytes_[i], keys_[i], 32);
    }
    initialized_ = true;
}

int64_t PRF::EvalMod64(KeyType type, int32_t index, int64_t p) const {
    if (!initialized_) {
        throw std::runtime_error("PRF not initialized: call KeyGen() first");
    }
    unsigned char msg[4];
    std::memcpy(msg, &index, 4);
    unsigned char digest[CryptoPP::HMAC<CryptoPP::SHA256>::DIGESTSIZE] = {0};
    CryptoPP::HMAC<CryptoPP::SHA256> hmac(key_bytes_[static_cast<int>(type)], 32);
    hmac.Update(msg, 4);
    hmac.Final(digest);
    // Interpret first 8 bytes of digest as uint64_t (little-endian).
    uint64_t x = 0;
    for (int i = 0; i < 8; ++i) {
        x |= static_cast<uint64_t>(digest[i]) << (8 * i);
    }
    // Reduce modulo p. Bias is ~p/2^64; negligible for p << 2^64.
    return static_cast<int64_t>(x % static_cast<uint64_t>(p));
}

NTL::ZZ PRF::Eval(KeyType type, int32_t index) const {
    if (!initialized_) {
        throw std::runtime_error("PRF not initialized: call KeyGen() first");
    }
    return HmacEval(keys_[static_cast<int>(type)], index);
}

NTL::ZZ PRF::EvalMod(KeyType type, int32_t index, const NTL::ZZ& p) const {
    return Eval(type, index) % p;
}

const NTL::ZZ& PRF::GetKey(KeyType type) const {
    return keys_[static_cast<int>(type)];
}

NTL::ZZ PRF::HmacEval(const NTL::ZZ& key, int32_t index) const {
    // Serialize the key to bytes
    unsigned char k[KEY_BYTES] = {0};
    NTL::BytesFromZZ(k, key, KEY_BYTES);

    // Serialize the index to bytes (little-endian)
    unsigned char m[MSG_BYTES] = {0};
    memcpy(m, &index, MSG_BYTES);

    // HMAC-SHA256
    unsigned char digest[CryptoPP::HMAC<CryptoPP::SHA256>::DIGESTSIZE] = {0};
    CryptoPP::HMAC<CryptoPP::SHA256> hmac(k, sizeof(k));
    hmac.Update(m, sizeof(m));
    hmac.Final(digest);

    // Convert digest to NTL::ZZ
    NTL::ZZ result;
    NTL::ZZFromBytes(result, digest, sizeof(digest));
    return result;
}

} // namespace vdpsu

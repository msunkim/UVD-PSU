// Test whether OpenFHE BFV supports a large plaintext modulus
// with N=2^13 for SIMD batching

#include <iostream>
#include <iomanip>
#include <vector>
#include "openfhe.h"

using namespace lbcrypto;

// Try different plaintext moduli and report which work
int main() {
    std::cout << "Testing large plaintext moduli with N=8192, depth=2" << std::endl;
    std::cout << "====================================================" << std::endl;

    // Candidate primes, each prime and = 1 mod 16384 for N=8192
    // These are primes of increasing size
    std::vector<uint64_t> candidates = {
        65537ULL,               // 17 bits (current)
        786433ULL,              // 20 bits
        1073872897ULL,          // 30 bits
        1099511922689ULL,       // 40 bits
        1125899906990081ULL,    // 50 bits
        576460752303439873ULL,  // 59 bits
    };

    for (uint64_t p : candidates) {
        std::cout << "\n--- Trying P = " << p
                  << " (log2 = " << std::fixed << std::setprecision(1)
                  << std::log2((double)p) << " bits) ---" << std::endl;

        // Check P = 1 mod 16384
        if ((p - 1) % 16384 != 0) {
            std::cout << "  SKIP: P != 1 mod 16384" << std::endl;
            continue;
        }

        try {
            // Auto-select N based on P size
            uint32_t N = (p > (1ULL << 40)) ? 16384 : 8192;
            CCParams<CryptoContextBFVRNS> params;
            params.SetRingDim(N);
            params.SetMultiplicativeDepth(2);
            params.SetPlaintextModulus(p);
            params.SetSecurityLevel(HEStd_128_classic);

            auto cc = GenCryptoContext(params);
            cc->Enable(PKE);
            cc->Enable(KEYSWITCH);
            cc->Enable(LEVELEDSHE);

            auto actual_N = cc->GetCryptoParameters()->GetElementParams()->GetRingDimension();
            auto batch = cc->GetEncodingParams()->GetBatchSize();
            std::cout << "  Ring dim: " << actual_N
                      << ", slots: " << batch << std::endl;

            // Test basic operation
            auto kp = cc->KeyGen();
            cc->EvalMultKeyGen(kp.secretKey);

            std::vector<int64_t> a(batch, 12345);
            std::vector<int64_t> b(batch, 67890);
            auto pt_a = cc->MakePackedPlaintext(a);
            auto pt_b = cc->MakePackedPlaintext(b);
            auto ct_a = cc->Encrypt(kp.publicKey, pt_a);
            auto ct_b = cc->Encrypt(kp.publicKey, pt_b);

            auto ct_mul = cc->EvalMult(ct_a, ct_b);
            auto ct_mul2 = cc->EvalMult(ct_mul, pt_b);  // depth 2

            Plaintext result;
            cc->Decrypt(kp.secretKey, ct_mul2, &result);
            result->SetLength(1);
            auto vals = result->GetPackedValue();

            // Expected: 12345 * 67890 * 67890 = 12345 * 67890^2
            __int128 expected = (__int128)12345 * 67890 * 67890;
            expected = expected % p;

            std::cout << "  Expected: " << (int64_t)expected
                      << ", Got: " << vals[0] << std::endl;

            if ((int64_t)expected == vals[0]) {
                std::cout << "  RESULT: OK" << std::endl;
            } else {
                std::cout << "  RESULT: WRONG VALUE" << std::endl;
            }
        } catch (const std::exception& e) {
            std::cout << "  RESULT: EXCEPTION: " << e.what() << std::endl;
        }
    }

    return 0;
}

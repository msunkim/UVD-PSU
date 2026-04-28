// Parameter comparison test: N=2^13 vs N=2^15
// Tests raw BFV operation costs at n=2^14 (m=2^15, representative scale)

#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include "openfhe.h"

using namespace lbcrypto;

struct ParamResult {
    size_t ring_dim;
    size_t slots;
    double encrypt_ms;
    double ct_ct_mul_ms;
    double pt_ct_mul_ms;
    double add_ms;
    double decrypt_ms;
    size_t ct_size_bytes;
};

ParamResult BenchmarkParams(uint32_t ring_dim, uint64_t plain_mod, uint32_t depth) {
    ParamResult r;
    r.ring_dim = ring_dim;

    CCParams<CryptoContextBFVRNS> params;
    params.SetRingDim(ring_dim);
    params.SetMultiplicativeDepth(depth);
    params.SetPlaintextModulus(plain_mod);
    params.SetSecurityLevel(HEStd_128_classic);

    auto cc = GenCryptoContext(params);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    auto actual_dim = cc->GetCryptoParameters()->GetElementParams()->GetRingDimension();
    auto batch = cc->GetEncodingParams()->GetBatchSize();
    r.slots = batch;

    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);

    // Prepare test vectors
    std::vector<int64_t> vec_a(batch, 1);
    std::vector<int64_t> vec_b(batch, 2);

    auto pt_a = cc->MakePackedPlaintext(vec_a);
    auto pt_b = cc->MakePackedPlaintext(vec_b);

    // Warmup
    auto ct_a = cc->Encrypt(kp.publicKey, pt_a);
    auto ct_b = cc->Encrypt(kp.publicKey, pt_b);

    using clk = std::chrono::high_resolution_clock;

    // Benchmark Encrypt
    {
        const int iters = 5;
        auto t0 = clk::now();
        for (int i = 0; i < iters; ++i) {
            auto ct = cc->Encrypt(kp.publicKey, pt_a);
        }
        auto t1 = clk::now();
        r.encrypt_ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / iters;
    }

    // Benchmark CT-CT multiply
    {
        const int iters = 5;
        auto t0 = clk::now();
        for (int i = 0; i < iters; ++i) {
            auto ct = cc->EvalMult(ct_a, ct_b);
        }
        auto t1 = clk::now();
        r.ct_ct_mul_ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / iters;
    }

    // Benchmark PT-CT multiply
    {
        const int iters = 5;
        auto t0 = clk::now();
        for (int i = 0; i < iters; ++i) {
            auto ct = cc->EvalMult(ct_a, pt_b);
        }
        auto t1 = clk::now();
        r.pt_ct_mul_ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / iters;
    }

    // Benchmark Add
    {
        const int iters = 20;
        auto t0 = clk::now();
        for (int i = 0; i < iters; ++i) {
            auto ct = cc->EvalAdd(ct_a, ct_b);
        }
        auto t1 = clk::now();
        r.add_ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / iters;
    }

    // Benchmark Decrypt
    {
        const int iters = 10;
        auto t0 = clk::now();
        for (int i = 0; i < iters; ++i) {
            Plaintext pt;
            cc->Decrypt(kp.secretKey, ct_a, &pt);
        }
        auto t1 = clk::now();
        r.decrypt_ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / iters;
    }

    // Estimate ciphertext size (serialized form)
    // ct has 2-3 ring elements, each of size ring_dim * log2(Q)/8 bytes
    // Rough estimate: use the ring dim and assume ~180 bits Q for depth 2
    r.ct_size_bytes = 2 * actual_dim * 23;  // ~23 bytes per coefficient (~180 bits)

    return r;
}

int main() {
    std::cout << "BFV Parameter Comparison: N=2^13 vs N=2^15" << std::endl;
    std::cout << "===========================================" << std::endl;

    // Use a prime that supports batching with both N values
    // For N=2^13: need P = 1 mod 16384
    // For N=2^15: need P = 1 mod 65536
    // 65537 = 1 mod 65536 (both should work)
    uint64_t plain_mod = 65537;
    uint32_t depth = 2;

    std::cout << "Plaintext modulus: " << plain_mod << std::endl;
    std::cout << "Multiplicative depth: " << depth << std::endl;
    std::cout << std::endl;

    std::cout << "Testing N = 2^13 = 8192..." << std::endl;
    auto r13 = BenchmarkParams(8192, plain_mod, depth);
    std::cout << "  Actual ring dim: " << r13.ring_dim
              << ", slots: " << r13.slots << std::endl;
    std::cout << "  Encrypt:  " << std::fixed << std::setprecision(3) << r13.encrypt_ms << " ms" << std::endl;
    std::cout << "  CT-CT x:  " << r13.ct_ct_mul_ms << " ms" << std::endl;
    std::cout << "  PT-CT x:  " << r13.pt_ct_mul_ms << " ms" << std::endl;
    std::cout << "  Add:      " << r13.add_ms << " ms" << std::endl;
    std::cout << "  Decrypt:  " << r13.decrypt_ms << " ms" << std::endl;
    std::cout << "  CT size:  " << r13.ct_size_bytes << " bytes" << std::endl;
    std::cout << std::endl;

    std::cout << "Testing N = 2^15 = 32768..." << std::endl;
    auto r15 = BenchmarkParams(32768, plain_mod, depth);
    std::cout << "  Actual ring dim: " << r15.ring_dim
              << ", slots: " << r15.slots << std::endl;
    std::cout << "  Encrypt:  " << std::fixed << std::setprecision(3) << r15.encrypt_ms << " ms" << std::endl;
    std::cout << "  CT-CT x:  " << r15.ct_ct_mul_ms << " ms" << std::endl;
    std::cout << "  PT-CT x:  " << r15.pt_ct_mul_ms << " ms" << std::endl;
    std::cout << "  Add:      " << r15.add_ms << " ms" << std::endl;
    std::cout << "  Decrypt:  " << r15.decrypt_ms << " ms" << std::endl;
    std::cout << "  CT size:  " << r15.ct_size_bytes << " bytes" << std::endl;
    std::cout << std::endl;

    // Estimate for n=2^20 (m=2^21)
    size_t m = 1 << 21;
    size_t cts_13 = (m + r13.slots - 1) / r13.slots;
    size_t cts_15 = (m + r15.slots - 1) / r15.slots;

    std::cout << "=== Estimated cost for n = 2^20 (m = 2^21) ===" << std::endl;
    std::cout << std::endl;
    std::cout << std::setw(20) << "" << std::setw(15) << "N=2^13" << std::setw(15) << "N=2^15" << std::endl;
    std::cout << std::setw(20) << "CTs per vector" << std::setw(15) << cts_13 << std::setw(15) << cts_15 << std::endl;

    // Compute phase: 3 CT-CT mul + 2 PT-CT mul + 2 add (on each of cts_X ciphertexts)
    double compute_13 = cts_13 * (3 * r13.ct_ct_mul_ms + 2 * r13.pt_ct_mul_ms + 2 * r13.add_ms);
    double compute_15 = cts_15 * (3 * r15.ct_ct_mul_ms + 2 * r15.pt_ct_mul_ms + 2 * r15.add_ms);
    std::cout << std::setw(20) << "Compute (sec)"
              << std::setw(15) << std::fixed << std::setprecision(1) << compute_13 / 1000.0
              << std::setw(15) << compute_15 / 1000.0 << std::endl;

    // Delegate phase: 5 vectors x cts encryptions
    double delegate_13 = 5 * cts_13 * r13.encrypt_ms;
    double delegate_15 = 5 * cts_15 * r15.encrypt_ms;
    std::cout << std::setw(20) << "Delegate (sec)"
              << std::setw(15) << delegate_13 / 1000.0
              << std::setw(15) << delegate_15 / 1000.0 << std::endl;

    // Decode phase: 4 vectors x cts decryptions
    double decode_13 = 4 * cts_13 * r13.decrypt_ms;
    double decode_15 = 4 * cts_15 * r15.decrypt_ms;
    std::cout << std::setw(20) << "Decode (sec)"
              << std::setw(15) << decode_13 / 1000.0
              << std::setw(15) << decode_15 / 1000.0 << std::endl;

    // Total
    double total_13 = compute_13 + delegate_13 + decode_13;
    double total_15 = compute_15 + delegate_15 + decode_15;
    std::cout << std::setw(20) << "Total (sec)"
              << std::setw(15) << total_13 / 1000.0
              << std::setw(15) << total_15 / 1000.0 << std::endl;

    // Communication
    size_t comm_13 = (5 + 4) * cts_13 * r13.ct_size_bytes;  // 5 up + 4 down
    size_t comm_15 = (5 + 4) * cts_15 * r15.ct_size_bytes;
    std::cout << std::setw(20) << "Comm (MB)"
              << std::setw(15) << comm_13 / (1024.0 * 1024.0)
              << std::setw(15) << comm_15 / (1024.0 * 1024.0) << std::endl;

    std::cout << std::endl;
    std::cout << "=== Winner ===" << std::endl;
    if (total_13 < total_15) {
        std::cout << "N = 2^13 is faster by "
                  << std::setprecision(1) << (100.0 * (total_15 - total_13) / total_15)
                  << "%" << std::endl;
    } else {
        std::cout << "N = 2^15 is faster by "
                  << std::setprecision(1) << (100.0 * (total_13 - total_15) / total_13)
                  << "%" << std::endl;
    }

    return 0;
}

#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>

#include "params.h"
#include "prf.h"
#include "mht.h"
#include "phf.h"
#include "timer.h"
#include "fhe_utils.h"
#include "client.h"
#include "cloud.h"

#include <set>
#include <algorithm>

using namespace vdpsu;

// ============================================================
// Test PRF
// ============================================================
void TestPRF() {
    std::cout << "=== TestPRF ===" << std::endl;

    PRF prf;
    prf.KeyGen(NTL::to_ZZ(12345L));

    // Determinism: same key + index = same output
    auto v1 = prf.Eval(KeyType::K1, 42);
    auto v2 = prf.Eval(KeyType::K1, 42);
    assert(v1 == v2);

    // Different indices produce different outputs
    auto v3 = prf.Eval(KeyType::K1, 43);
    assert(v1 != v3);

    // Different key types produce different outputs
    auto v4 = prf.Eval(KeyType::K2, 42);
    assert(v1 != v4);

    std::cout << "  PRF determinism: PASS" << std::endl;
    std::cout << "  PRF distinct outputs: PASS" << std::endl;
}

// ============================================================
// Test MHT
// ============================================================
void TestMHT() {
    std::cout << "=== TestMHT ===" << std::endl;

    MHT mht(10);
    assert(mht.Size() == 10);
    assert(mht.Query(3) == 0);

    mht.Insert(3);
    assert(mht.Query(3) == 1);

    mht.Remove(3);
    assert(mht.Query(3) == 0);

    // Hadamard product
    MHT a(4), b(4);
    a.Insert(0); a.Insert(2);  // [1, 0, 1, 0]
    std::vector<int64_t> vec = {5, 7, 3, 9};
    auto prod = a.HadamardProduct(vec);
    assert(prod[0] == 5 && prod[1] == 0 && prod[2] == 3 && prod[3] == 0);

    std::cout << "  MHT basic ops: PASS" << std::endl;
    std::cout << "  MHT Hadamard: PASS" << std::endl;
}

// ============================================================
// Test PHF
// ============================================================
void TestPHF() {
    std::cout << "=== TestPHF ===" << std::endl;

    std::vector<NTL::ZZ> elements;
    for (int i = 0; i < 5; ++i) {
        elements.push_back(NTL::to_ZZ(static_cast<long>(100 + i)));
    }

    PHF phf;
    phf.Build(elements, 10);

    // Injectivity: all indices should be distinct
    std::set<size_t> indices;
    for (const auto& e : elements) {
        size_t idx = phf.Eval(e);
        assert(idx < 10);
        indices.insert(idx);
    }
    assert(indices.size() == 5);  // all distinct

    // Inverse
    for (const auto& e : elements) {
        size_t idx = phf.Eval(e);
        assert(phf.InverseEval(idx) == e);
    }

    std::cout << "  PHF injectivity: PASS" << std::endl;
    std::cout << "  PHF inverse: PASS" << std::endl;
}

// ============================================================
// Test FHE Utils (BFV SIMD)
// ============================================================
void TestFHEUtils() {
    std::cout << "=== TestFHEUtils ===" << std::endl;

    FHEUtils fhe;
    fhe.Init();  // default plaintext modulus = 65537

    auto kp = fhe.KeyGen();

    // Test 1: encrypt-decrypt roundtrip (small vector)
    {
        std::vector<int64_t> vec = {1, 0, 1, 0, 1, 1, 0, 0};
        auto ct = fhe.Encrypt(kp, vec);
        auto dec = fhe.Decrypt(kp, ct, vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            assert(dec[i] == vec[i]);
        }
        std::cout << "  FHE roundtrip (small): PASS" << std::endl;
    }

    // Test 2: encrypt-decrypt roundtrip (larger vector, multiple ciphertexts)
    {
        size_t m = 5000;  // > 4096, requires 2 ciphertexts
        std::vector<int64_t> vec(m);
        for (size_t i = 0; i < m; ++i) vec[i] = i % 3;

        auto ct = fhe.Encrypt(kp, vec);
        std::cout << "  num ciphertexts for m=5000: " << ct.size() << std::endl;

        auto dec = fhe.Decrypt(kp, ct, m);
        for (size_t i = 0; i < m; ++i) {
            assert(dec[i] == vec[i]);
        }
        std::cout << "  FHE roundtrip (multi-ct): PASS" << std::endl;
    }

    // Test 3: homomorphic addition
    {
        std::vector<int64_t> a = {1, 0, 1, 0};
        std::vector<int64_t> b = {0, 1, 1, 0};
        auto ct_a = fhe.Encrypt(kp, a);
        auto ct_b = fhe.Encrypt(kp, b);
        auto ct_sum = fhe.HAdd(ct_a, ct_b);
        auto dec = fhe.Decrypt(kp, ct_sum, 4);
        assert(dec[0] == 1 && dec[1] == 1 && dec[2] == 2 && dec[3] == 0);
        std::cout << "  FHE HAdd: PASS" << std::endl;
    }

    // Test 4: ciphertext-ciphertext multiplication (depth-1)
    {
        std::vector<int64_t> a = {1, 0, 1, 0, 1};
        std::vector<int64_t> b = {1, 1, 0, 1, 1};
        auto ct_a = fhe.Encrypt(kp, a);
        auto ct_b = fhe.Encrypt(kp, b);
        auto ct_prod = fhe.HMul(ct_a, ct_b);
        auto dec = fhe.Decrypt(kp, ct_prod, 5);
        // Expected: [1*1, 0*1, 1*0, 0*1, 1*1] = [1, 0, 0, 0, 1]
        assert(dec[0] == 1 && dec[1] == 0 && dec[2] == 0 && dec[3] == 0 && dec[4] == 1);
        std::cout << "  FHE HMul (ct-ct): PASS" << std::endl;
    }

    // Test 5: plaintext-ciphertext multiplication
    {
        std::vector<int64_t> pt = {3, 0, 5, 7};
        std::vector<int64_t> ct_vals = {1, 1, 0, 1};
        auto ct = fhe.Encrypt(kp, ct_vals);
        auto ct_prod = fhe.PTMul(pt, ct);
        auto dec = fhe.Decrypt(kp, ct_prod, 4);
        // Expected: [3*1, 0*1, 5*0, 7*1] = [3, 0, 0, 7]
        assert(dec[0] == 3 && dec[1] == 0 && dec[2] == 0 && dec[3] == 7);
        std::cout << "  FHE PTMul: PASS" << std::endl;
    }

    // Test 6: simulate the Compute phase formula
    // delta = tilde_mht_B + (-1)*mht_A
    // mht = tilde_mht_B * delta
    {
        // Suppose A = {x1, x3}, B = {x2, x3} in positions [0..4]
        // mht_A = [1, 0, 1, 0, 0], tilde_mht_A = [1, 0, 1, 0, 0]
        // mht_B = [0, 1, 1, 0, 0], tilde_mht_B = [0, 1, 1, 0, 0]
        std::vector<int64_t> neg_mht_A = {-1, 0, -1, 0, 0};  // (-1)*mht_A
        std::vector<int64_t> tilde_mht_B = {0, 1, 1, 0, 0};

        auto ct_neg_mht_A = fhe.Encrypt(kp, neg_mht_A);
        auto ct_tilde_mht_B = fhe.Encrypt(kp, tilde_mht_B);

        // delta = tilde_mht_B + (-1)*mht_A
        auto ct_delta = fhe.HAdd(ct_tilde_mht_B, ct_neg_mht_A);
        auto delta_dec = fhe.Decrypt(kp, ct_delta, 5);
        // Expected: [0+(-1), 1+0, 1+(-1), 0+0, 0+0] = [-1, 1, 0, 0, 0]
        // BFV may represent -1 as either -1 or p-1 depending on OpenFHE's decoding
        std::cout << "  delta[0] = " << delta_dec[0]
                  << " (expected -1 or 65536)" << std::endl;
        assert(delta_dec[0] == -1 || delta_dec[0] == 65536);
        assert(delta_dec[1] == 1);
        assert(delta_dec[2] == 0);

        // mht = tilde_mht_B * delta (depth-1 ct-ct multiply)
        auto ct_mht = fhe.HMul(ct_tilde_mht_B, ct_delta);
        auto mht_dec = fhe.Decrypt(kp, ct_mht, 5);
        // Expected: [0*(-1), 1*1, 1*0, 0*0, 0*0] = [0, 1, 0, 0, 0]
        assert(mht_dec[1] == 1);
        assert(mht_dec[2] == 0);
        assert(mht_dec[3] == 0);
        std::cout << "  FHE Compute simulation: PASS" << std::endl;
    }
}

// ============================================================
// Test End-to-End VD-PSU Protocol
// ============================================================
void TestProtocolE2E() {
    std::cout << "=== TestProtocolE2E ===" << std::endl;

    // X_A = {10, 20, 30, 40}, X_B = {20, 30, 40, 50}
    // Expected union = {10, 20, 30, 40, 50}
    std::vector<int64_t> X_A = {10, 20, 30, 40};
    std::vector<int64_t> X_B = {20, 30, 40, 50};
    std::set<int64_t> expected_union = {10, 20, 30, 40, 50};

    size_t n = 4;
    ProtocolConfig cfg(n);
    std::cout << "  n=" << cfg.n << ", m=" << cfg.m
              << ", num_cts=" << cfg.num_cts << std::endl;

    // Shared FHE context
    FHEUtils fhe;
    fhe.Init();  // plaintext modulus = 65537

    // --- Setup: Build the shared PHF (public parameter) ---
    // Both clients first compute their witnesses, then we build a PHF
    // over the union of all elements: X_A + X_B + wit_A + wit_B.
    // In a real deployment, the PHF would be built from a universal domain.
    // Here we simulate the setup phase.

    // Temporarily create clients to get their witnesses
    // (In practice, the PHF is a public parameter agreed upon before encoding.)
    // We need a two-pass approach: first get witnesses, then build PHF, then encode.
    PHF shared_phf;

    // Create clients with a placeholder PHF first to get witnesses
    // Actually, we construct the clients, call ComputeWitness, build PHF, then Encode.
    ClientA clientA(cfg, fhe, shared_phf);
    ClientB clientB(cfg, fhe, shared_phf);

    NTL::ZZ wit_A = clientA.ComputeWitness();
    NTL::ZZ wit_B = clientB.ComputeWitness();

    // Build shared PHF from all elements: X_A union X_B union {wit_A, wit_B}
    std::set<int64_t> all_elems_set;
    for (auto x : X_A) all_elems_set.insert(x);
    for (auto x : X_B) all_elems_set.insert(x);
    all_elems_set.insert(NTL::conv<long>(wit_A));
    all_elems_set.insert(NTL::conv<long>(wit_B));

    std::vector<NTL::ZZ> all_elems;
    for (auto x : all_elems_set) {
        all_elems.push_back(NTL::to_ZZ(static_cast<long>(x)));
    }
    shared_phf.Build(all_elems, cfg.m);

    std::cout << "  Shared PHF built over " << all_elems.size()
              << " elements into m=" << cfg.m << " slots" << std::endl;

    // Create cloud
    Cloud cloud(cfg, fhe);

    Timer timer;

    // Phase 1: Encode
    timer.Start();
    auto bar_e_A = clientA.Encode(X_A);
    auto bar_e_B = clientB.Encode(X_B);
    double encode_ms = timer.StopMs();
    std::cout << "  Encode: " << encode_ms << " ms" << std::endl;

    assert(bar_e_A.size() == cfg.m);
    assert(bar_e_B.size() == cfg.m);
    std::cout << "  Encode output sizes: PASS" << std::endl;

    // Phase 2: Delegate
    timer.Start();
    auto token_A_cloud = clientA.Delegate();
    auto token_B_cloud = clientB.Delegate(clientA.GetKeyPair());
    auto token_B_to_A = clientB.GetTokenForA();
    double delegate_ms = timer.StopMs();
    std::cout << "  Delegate: " << delegate_ms << " ms" << std::endl;

    // Phase 3: Cloud Compute
    timer.Start();
    auto cloud_response = cloud.Compute(token_A_cloud, token_B_cloud,
                                         bar_e_A, bar_e_B);
    double compute_ms = timer.StopMs();
    std::cout << "  Compute: " << compute_ms << " ms" << std::endl;

    // Phase 4: Decode + Verify
    timer.Start();
    auto union_result = clientA.Decode(cloud_response, token_B_to_A);
    double decode_ms = timer.StopMs();
    std::cout << "  Decode: " << decode_ms << " ms" << std::endl;

    // Check result
    assert(!union_result.empty() && "Decode returned empty (verification failure)");

    std::set<int64_t> result_set(union_result.begin(), union_result.end());

    std::cout << "  Union result: {";
    for (auto it = result_set.begin(); it != result_set.end(); ++it) {
        if (it != result_set.begin()) std::cout << ", ";
        std::cout << *it;
    }
    std::cout << "}" << std::endl;

    std::cout << "  Expected:     {";
    for (auto it = expected_union.begin(); it != expected_union.end(); ++it) {
        if (it != expected_union.begin()) std::cout << ", ";
        std::cout << *it;
    }
    std::cout << "}" << std::endl;

    assert(result_set == expected_union && "Union does not match expected!");
    std::cout << "  End-to-end protocol: PASS" << std::endl;

    double total_ms = encode_ms + delegate_ms + compute_ms + decode_ms;
    std::cout << "  Total time: " << total_ms << " ms" << std::endl;
}

// ============================================================
// Main
// ============================================================
int main() {
    std::cout << "VD-PSU Unit Tests" << std::endl;
    std::cout << "=================" << std::endl;

    TestPRF();
    TestMHT();
    TestPHF();
    TestFHEUtils();
    TestProtocolE2E();

    std::cout << std::endl;
    std::cout << "All tests PASSED." << std::endl;
    return 0;
}

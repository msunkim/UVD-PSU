#include "client.h"

#include <algorithm>
#include <atomic>
#include <iostream>
#include <random>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_thread_num() { return 0; }
#endif

namespace vdpsu {

// ============================================================
// ClientA helpers
// ============================================================
int64_t ClientA::ModP(int64_t a) {
    int64_t r = a % P;
    if (r < 0) r += P;
    return r;
}

// Safe modular multiplication using __int128 to avoid int64_t overflow
static inline int64_t MulModP_impl(int64_t a, int64_t b, int64_t P) {
    __int128 prod = static_cast<__int128>(a) * static_cast<__int128>(b);
    int64_t r = static_cast<int64_t>(prod % P);
    if (r < 0) r += P;
    return r;
}

int64_t ClientA::MulModP(int64_t a, int64_t b) {
    return MulModP_impl(a, b, P);
}

int64_t ClientA::ModInverse(int64_t a) {
    a = ModP(a);
    // Extended Euclidean algorithm using __int128 to guard against overflow
    // in the `q * y` update when P approaches the int64 range.
    __int128 g = P, x = 0, y = 1;
    __int128 a0 = a;
    while (a0 != 0) {
        __int128 q = g / a0;
        __int128 tmp = g - q * a0;
        g = a0; a0 = tmp;
        tmp = x - q * y;
        x = y; y = tmp;
    }
    if (g != 1) {
        throw std::runtime_error("ModInverse: no inverse exists");
    }
    int64_t r = static_cast<int64_t>(x % P);
    if (r < 0) r += P;
    return r;
}

// ============================================================
// ClientB helpers
// ============================================================
int64_t ClientB::MulModP(int64_t a, int64_t b) {
    return MulModP_impl(a, b, P);
}

int64_t ClientB::ModP(int64_t a) {
    int64_t r = a % P;
    if (r < 0) r += P;
    return r;
}

// ============================================================
// ClientA
// ============================================================
ClientA::ClientA(const ProtocolConfig& cfg, FHEUtils& fhe, const PHF& phf)
    : cfg_(cfg), fhe_(fhe), phf_(phf), tilde_mht_(cfg.m),
      bf_A_(cfg.n + 1, BF_NUM_HASHES),
      mask_bf_A_(cfg.n + 1, BF_NUM_HASHES),
      xmask_bf_A_(cfg.n + 1, BF_NUM_HASHES) {
    prf_.KeyGen();
}

NTL::ZZ ClientA::ComputeWitness() {
    wit_A_ = prf_.EvalMod(KeyType::K3, 0, NTL::to_ZZ(static_cast<long>(P)));
    wit_A_val_ = NTL::conv<long>(wit_A_);
    return wit_A_;
}

std::vector<int64_t> ClientA::Encode(const std::vector<int64_t>& X_A) {
    size_t n = X_A.size();
    size_t m = cfg_.m;

    // Build tilde_mht: insert all elements of X_A and wit_A using shared PHF.
    // Each PHF.Eval is O(1) hashmap lookup; MHT.Insert is a single byte write.
    // Sequential is fine here — this loop is << 5% of Encode time.
    tilde_mht_ = MHT(m);
    for (size_t i = 0; i < n; ++i) {
        size_t idx = phf_.Eval(X_A[i]);
        tilde_mht_.Insert(idx);
    }
    tilde_mht_.Insert(phf_.Eval(wit_A_));

    // Main encoding loop — parallelized. Writes to bar_e[j] are safe because
    // PHF is injective (each iteration writes a distinct index), and BF inserts
    // use atomic fetch-or internally.
    std::vector<int64_t> bar_e(m, 0);

    #pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        size_t j = phf_.Eval(X_A[i]);
        int64_t r_j = prf_.EvalMod64(KeyType::K1, static_cast<int32_t>(j), P);
        int64_t s_j = prf_.EvalMod64(KeyType::K2, static_cast<int32_t>(j), P);

        int64_t bar_x = ModP(MulModP(r_j, X_A[i]) + s_j);
        bar_e[j] = bar_x;

        bf_A_.Insert(X_A[i]);
        mask_bf_A_.Insert(bar_x);

        int64_t xmask_x = ModP(MulModP(s_j, X_A[i]) + r_j);
        xmask_bf_A_.Insert(xmask_x);
    }

    // Witness encoding (sequential — single element).
    {
        size_t j = phf_.Eval(wit_A_);
        int64_t r_j = prf_.EvalMod64(KeyType::K1, static_cast<int32_t>(j), P);
        int64_t s_j = prf_.EvalMod64(KeyType::K2, static_cast<int32_t>(j), P);

        int64_t bar_wit = ModP(MulModP(r_j, wit_A_val_) + s_j);
        bar_e[j] = bar_wit;
        mask_bf_A_.Insert(bar_wit);

        int64_t xmask_wit = ModP(MulModP(s_j, wit_A_val_) + r_j);
        xmask_bf_A_.Insert(xmask_wit);
    }

    // Fill wildcard (unoccupied) positions with random dummy values.
    // Parallelized with per-thread RNG so no shared state.
    #pragma omp parallel
    {
        std::mt19937_64 rng(static_cast<uint64_t>(std::random_device{}()) +
                            static_cast<uint64_t>(omp_get_thread_num()));
        std::uniform_int_distribution<int64_t> dist(1, P - 1);
        #pragma omp for nowait
        for (size_t j = 0; j < m; ++j) {
            if (tilde_mht_.Query(j) == 0) {
                bar_e[j] = dist(rng);
            }
        }
    }

    return bar_e;
}

TokenA_Cloud ClientA::Delegate() {
    size_t m = cfg_.m;

    // Generate FHE key pair
    kp_ = fhe_.KeyGen();

    // tilde_mht_ is mht_A (with witness).
    // Compute tilde_mht_A = mht_A with wit_A removed.
    MHT tilde_mht_A(m);
    auto& mht_A_vec = tilde_mht_.GetVector();
    for (size_t j = 0; j < m; ++j) {
        if (mht_A_vec[j] == 1) tilde_mht_A.Insert(j);
    }
    size_t wit_idx = phf_.Eval(wit_A_);
    tilde_mht_A.Remove(wit_idx);

    // Encrypt mht_A (with witness) and (-1)*tilde_mht_A (without witness)
    MHT neg_tilde_mht_A = tilde_mht_A.Negate();

    TokenA_Cloud token;
    token.ct_mht_A = fhe_.Encrypt(kp_, tilde_mht_.GetVector());
    token.ct_neg_tilde_mht_A = fhe_.Encrypt(kp_, neg_tilde_mht_A.GetVector());
    return token;
}

std::vector<int64_t> ClientA::Decode(const CloudResponse& response,
                                     const TokenB_ToA& token_B) {
    size_t m = cfg_.m;

    // Decrypt all four ciphertext vectors
    auto bar_y     = fhe_.Decrypt(kp_, response.ct_bar_y, m);
    auto vecd_vec  = fhe_.Decrypt(kp_, response.ct_vecd, m);
    auto xmask_r_B = fhe_.Decrypt(kp_, response.ct_xmask_r_B, m);
    auto xmask_s_B = fhe_.Decrypt(kp_, response.ct_xmask_s_B, m);

    // Normalize decrypted values to [0, P-1] (parallel).
    #pragma omp parallel for
    for (size_t j = 0; j < m; ++j) {
        bar_y[j]     = ModP(bar_y[j]);
        vecd_vec[j]  = ModP(vecd_vec[j]);
        xmask_r_B[j] = ModP(xmask_r_B[j]);
        xmask_s_B[j] = ModP(xmask_s_B[j]);
    }

    // Compute union indicator mht_U = mht_A + vecd
    // tilde_mht_ holds mht_A (with witness). vecd is the difference indicator.
    // mht_U = mht_A + vecd (arithmetic addition, per Algorithm 4)
    // For honest execution: mht_U[j] in {0, 1}.
    // A malicious cloud producing vecd[j] > 1 would yield mht_U[j] > 1, detected below.
    auto mht_A_vec = tilde_mht_.GetVector();
    std::vector<int64_t> mht_U(m);
    #pragma omp parallel for
    for (size_t j = 0; j < m; ++j) {
        mht_U[j] = mht_A_vec[j] + vecd_vec[j];
    }

    // Consistency check: mht_U[j] must be 0 or 1, and mht_U[j]=1 iff bar_y[j]!=0
    int bad_mhtU = 0, bad_consistency = 0;
    int64_t first_bad_mhtU_val = 0;
    size_t first_bad_mhtU_j = 0, first_bad_cons_j = 0;
    for (size_t j = 0; j < m; ++j) {
        if (mht_U[j] < 0 || mht_U[j] > 1) {
            if (bad_mhtU == 0) { first_bad_mhtU_val = mht_U[j]; first_bad_mhtU_j = j; }
            bad_mhtU++;
        }
        bool has_elem = (mht_U[j] == 1);
        bool y_nonzero = (bar_y[j] != 0);
        if (has_elem != y_nonzero) {
            if (bad_consistency == 0) first_bad_cons_j = j;
            bad_consistency++;
        }
    }
    if (bad_mhtU > 0 || bad_consistency > 0) {
        std::cerr << "  [Decode FAIL m=" << m
                  << "] bad_mhtU=" << bad_mhtU
                  << " (first j=" << first_bad_mhtU_j << " val=" << first_bad_mhtU_val << ")"
                  << " bad_cons=" << bad_consistency
                  << " (first j=" << first_bad_cons_j << ")"
                  << std::endl;
        return {};
    }

    // Build tilde_bf_A = bf_A union {wit_A}
    BloomFilter tilde_bf_A = bf_A_;
    tilde_bf_A.Insert(wit_A_val_);

    // Recompute PRF vectors for A (parallel — each j is independent)
    std::vector<int64_t> r_A(m), s_A(m);
    #pragma omp parallel for
    for (size_t j = 0; j < m; ++j) {
        r_A[j] = prf_.EvalMod64(KeyType::K1, static_cast<int32_t>(j), P);
        s_A[j] = prf_.EvalMod64(KeyType::K2, static_cast<int32_t>(j), P);
    }

    // Recover the union in parallel. Each position j is independent; we write
    // into a dense per-slot array and compact to U in a sequential pass.
    // -1 marks an empty slot (no element at j) or an undecodable position.
    // A flag signals verification failure to all threads via atomic store.
    std::vector<int64_t> slot(m, -1);
    std::atomic<bool> verify_failed{false};

    #pragma omp parallel for
    for (size_t j = 0; j < m; ++j) {
        if (verify_failed.load(std::memory_order_relaxed)) continue;
        if (mht_U[j] == 0) continue;

        int64_t y_j = bar_y[j];
        int64_t x;

        if (mask_bf_A_.Query(y_j)) {
            // Element from A: bar_y_j = r_A[j]*x + s_A[j]
            int64_t num = ModP(y_j - s_A[j]);
            x = MulModP(num, ModInverse(r_A[j]));
            if (!tilde_bf_A.Query(x)) {
                verify_failed.store(true, std::memory_order_relaxed);
                continue;
            }
        } else {
            // Element from B\A: bar_y_j = r_B[j]*x + s_B[j]
            int64_t r_B_j = xmask_r_B[j];
            int64_t s_B_j = xmask_s_B[j];
            if (r_B_j == 0) {
                verify_failed.store(true, std::memory_order_relaxed);
                continue;
            }
            int64_t num = ModP(y_j - s_B_j);
            x = MulModP(num, ModInverse(r_B_j));
            // Verify: xmask_x = s_B[j]*x + r_B[j]
            int64_t xmask_x = ModP(MulModP(s_B_j, x) + r_B_j);
            if (!token_B.xmask_bf_B.Query(xmask_x)) {
                verify_failed.store(true, std::memory_order_relaxed);
                continue;
            }
        }
        slot[j] = x;
    }

    if (verify_failed.load()) return {};

    // Compact non-empty slots into U.
    std::vector<int64_t> U;
    U.reserve(m);
    for (size_t j = 0; j < m; ++j) {
        if (slot[j] != -1) U.push_back(slot[j]);
    }

    // Verify both witnesses are present
    int64_t wit_B_val = NTL::conv<long>(token_B.wit_B);
    bool found_A = false, found_B = false;
    for (auto v : U) {
        if (v == wit_A_val_) found_A = true;
        if (v == wit_B_val) found_B = true;
    }
    if (!found_A || !found_B) {
        return {};
    }

    // Remove witnesses
    std::vector<int64_t> result;
    for (auto v : U) {
        if (v != wit_A_val_ && v != wit_B_val) {
            result.push_back(v);
        }
    }
    return result;
}

// ============================================================
// ClientA::Update (Algorithm 5, party = A)
// ============================================================
bool ClientA::Update(char op, int64_t x, UpdateRequest& out_req) {
    // Flat-MHT prototype: bin index j = 0; pos = phf(x).
    size_t i = phf_.Eval(x);

    int64_t r = prf_.EvalMod64(KeyType::K1, static_cast<int32_t>(i), P);
    int64_t s = prf_.EvalMod64(KeyType::K2, static_cast<int32_t>(i), P);
    int64_t mask_x  = ModP(MulModP(r, x) + s);
    int64_t xmask_x = ModP(MulModP(s, x) + r);

    if (op == '+') {
        if (tilde_mht_.Query(i) == 1) return false;          // duplicate
        tilde_mht_.Insert(i);
        bf_A_.Insert(x);
        mask_bf_A_.Insert(mask_x);
        xmask_bf_A_.Insert(xmask_x);
        out_req = {'A', mask_x, /*bin*/ 0, i};
        return true;
    } else if (op == '-') {
        if (tilde_mht_.Query(i) == 0) return false;          // non-existent
        tilde_mht_.Remove(i);
        // Standard BFs are not deletable; in deployment a counting/dynamic BF
        // (Wang et al., WHY+22) would replace these. Skipping here keeps the
        // benchmark simple and does not affect timing of the dominant work.
        thread_local std::mt19937_64 rng(static_cast<uint64_t>(std::random_device{}()));
        std::uniform_int_distribution<int64_t> dist(1, P - 1);
        out_req = {'A', dist(rng), /*bin*/ 0, i};
        return true;
    }
    return false;
}

bool ClientA::RefreshWitness(int64_t new_witness,
                              UpdateRequest& del_req,
                              UpdateRequest& ins_req) {
    if (!Update('-', wit_A_val_, del_req)) return false;
    wit_A_val_ = new_witness;
    wit_A_     = NTL::to_ZZ(static_cast<long>(new_witness));
    return Update('+', new_witness, ins_req);
}

// ============================================================
// ClientB
// ============================================================
ClientB::ClientB(const ProtocolConfig& cfg, FHEUtils& fhe, const PHF& phf)
    : cfg_(cfg), fhe_(fhe), phf_(phf), tilde_mht_(cfg.m),
      bf_B_(cfg.n + 1, BF_NUM_HASHES),
      mask_bf_B_(cfg.n + 1, BF_NUM_HASHES),
      xmask_bf_B_(cfg.n + 1, BF_NUM_HASHES) {
    prf_.KeyGen();
}

NTL::ZZ ClientB::ComputeWitness() {
    wit_B_ = prf_.EvalMod(KeyType::K3, 0, NTL::to_ZZ(static_cast<long>(P)));
    wit_B_val_ = NTL::conv<long>(wit_B_);
    return wit_B_;
}

std::vector<int64_t> ClientB::Encode(const std::vector<int64_t>& X_B) {
    size_t n = X_B.size();
    size_t m = cfg_.m;

    // Build tilde_mht using shared PHF
    tilde_mht_ = MHT(m);
    for (size_t i = 0; i < n; ++i) {
        size_t idx = phf_.Eval(X_B[i]);
        tilde_mht_.Insert(idx);
    }
    tilde_mht_.Insert(phf_.Eval(wit_B_));

    // Compute masked encodings (parallelized — BF.Insert uses atomic fetch-or).
    std::vector<int64_t> bar_e(m, 0);

    #pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        size_t j = phf_.Eval(X_B[i]);
        int64_t r_j = prf_.EvalMod64(KeyType::K1, static_cast<int32_t>(j), P);
        int64_t s_j = prf_.EvalMod64(KeyType::K2, static_cast<int32_t>(j), P);

        int64_t bar_x = ModP(MulModP(r_j, X_B[i]) + s_j);
        bar_e[j] = bar_x;

        bf_B_.Insert(X_B[i]);
        mask_bf_B_.Insert(bar_x);

        int64_t xmask_x = ModP(MulModP(s_j, X_B[i]) + r_j);
        xmask_bf_B_.Insert(xmask_x);
    }

    // Witness
    {
        size_t j = phf_.Eval(wit_B_);
        int64_t r_j = prf_.EvalMod64(KeyType::K1, static_cast<int32_t>(j), P);
        int64_t s_j = prf_.EvalMod64(KeyType::K2, static_cast<int32_t>(j), P);

        int64_t bar_wit = ModP(MulModP(r_j, wit_B_val_) + s_j);
        bar_e[j] = bar_wit;
        mask_bf_B_.Insert(bar_wit);

        int64_t xmask_wit = ModP(MulModP(s_j, wit_B_val_) + r_j);
        xmask_bf_B_.Insert(xmask_wit);
    }

    // Fill wildcard positions (parallel, per-thread RNG).
    #pragma omp parallel
    {
        std::mt19937_64 rng(static_cast<uint64_t>(std::random_device{}()) +
                            static_cast<uint64_t>(omp_get_thread_num()) + 1);
        std::uniform_int_distribution<int64_t> dist(1, P - 1);
        #pragma omp for nowait
        for (size_t j = 0; j < m; ++j) {
            if (tilde_mht_.Query(j) == 0) {
                bar_e[j] = dist(rng);
            }
        }
    }

    return bar_e;
}

TokenB_Cloud ClientB::Delegate(const KeyPair& pk_from_A) {
    size_t m = cfg_.m;

    // Build r_B and s_B vectors (parallel — each j is independent)
    std::vector<int64_t> r_B(m), s_B(m);
    #pragma omp parallel for
    for (size_t j = 0; j < m; ++j) {
        r_B[j] = prf_.EvalMod64(KeyType::K1, static_cast<int32_t>(j), P);
        s_B[j] = prf_.EvalMod64(KeyType::K2, static_cast<int32_t>(j), P);
    }

    // tilde_mht_ is mht_B (with witness)
    auto mht_B_vec = tilde_mht_.GetVector();

    // Encrypt mht_B, r_B, s_B under A's public key
    TokenB_Cloud token;
    token.ct_mht_B = fhe_.Encrypt(pk_from_A, mht_B_vec);
    token.ct_r_B   = fhe_.Encrypt(pk_from_A, r_B);
    token.ct_s_B   = fhe_.Encrypt(pk_from_A, s_B);
    return token;
}

TokenB_ToA ClientB::GetTokenForA() const {
    TokenB_ToA token;
    token.xmask_bf_B = xmask_bf_B_;
    token.wit_B = wit_B_;
    return token;
}

// ============================================================
// ClientB::Update (Algorithm 5, party = B)
// ============================================================
bool ClientB::Update(char op, int64_t x, UpdateRequest& out_req) {
    size_t i = phf_.Eval(x);

    int64_t r = prf_.EvalMod64(KeyType::K1, static_cast<int32_t>(i), P);
    int64_t s = prf_.EvalMod64(KeyType::K2, static_cast<int32_t>(i), P);
    int64_t mask_x  = ModP(MulModP(r, x) + s);
    int64_t xmask_x = ModP(MulModP(s, x) + r);

    if (op == '+') {
        if (tilde_mht_.Query(i) == 1) return false;
        tilde_mht_.Insert(i);
        bf_B_.Insert(x);
        mask_bf_B_.Insert(mask_x);
        xmask_bf_B_.Insert(xmask_x);
        out_req = {'B', mask_x, /*bin*/ 0, i};
        return true;
    } else if (op == '-') {
        if (tilde_mht_.Query(i) == 0) return false;
        tilde_mht_.Remove(i);
        thread_local std::mt19937_64 rng(static_cast<uint64_t>(std::random_device{}()) ^ 0xB);
        std::uniform_int_distribution<int64_t> dist(1, P - 1);
        out_req = {'B', dist(rng), /*bin*/ 0, i};
        return true;
    }
    return false;
}

bool ClientB::RefreshWitness(int64_t new_witness,
                              UpdateRequest& del_req,
                              UpdateRequest& ins_req) {
    if (!Update('-', wit_B_val_, del_req)) return false;
    wit_B_val_ = new_witness;
    wit_B_     = NTL::to_ZZ(static_cast<long>(new_witness));
    return Update('+', new_witness, ins_req);
}

} // namespace vdpsu

#ifndef VDPSU_CLIENT_H
#define VDPSU_CLIENT_H

#include "params.h"
#include "prf.h"
#include "mht.h"
#include "phf.h"
#include "fhe_utils.h"
#include "bloom_filter.h"

#include <vector>
#include <string>
#include <NTL/ZZ.h>

namespace vdpsu {

// Token from Client A to cloud
struct TokenA_Cloud {
    CTVector ct_mht_A;           // Enc(mht_A) - with witness
    CTVector ct_neg_tilde_mht_A; // Enc(-tilde_mht_A) - without witness, negated
};

// Token from Client B to cloud
struct TokenB_Cloud {
    CTVector ct_mht_B;  // Enc(mht_B) - with witness
    CTVector ct_r_B;    // Enc(r_B)
    CTVector ct_s_B;    // Enc(s_B)
};

// Token from Client B to Client A
struct TokenB_ToA {
    BloomFilter xmask_bf_B;  // cross-mask BF of B
    NTL::ZZ wit_B;           // witness of B
};

// Cloud response to Client A
struct CloudResponse {
    CTVector ct_bar_y;     // encrypted union encoding
    CTVector ct_vecd;      // encrypted difference indicator
    CTVector ct_xmask_r_B; // encrypted masked r_B
    CTVector ct_xmask_s_B; // encrypted masked s_B
};

// An update request that a client sends to the cloud (Algorithm 5).
// Fields: party identifier ('A' or 'B'), masked value, bin index, permuted slot.
// In the flat-MHT prototype, the bin index is always 0 and `pos` is the global
// position in the size-m masked vector.
struct UpdateRequest {
    char    party;
    int64_t mask_x;
    size_t  bin;
    size_t  pos;
};

// ============================================================
// ClientA
// ============================================================
class ClientA {
public:
    // cfg: protocol config, fhe: shared FHE utils, phf: shared PHF (public param)
    ClientA(const ProtocolConfig& cfg, FHEUtils& fhe, const PHF& phf);

    // Compute the witness value wit_A = PRF(K3, ctr) mod P.
    // Must be called before Encode so the PHF can be built with the witness.
    NTL::ZZ ComputeWitness();

    // Phase 1: Encode input set X_A.
    // Builds MHT, masked encodings (bar_e_A), 3 BFs.
    // ComputeWitness() must have been called first.
    // Returns bar_e_A (the permuted masked encoding vector of length m).
    std::vector<int64_t> Encode(const std::vector<int64_t>& X_A);

    // Phase 2: Delegate.
    // Generates BFV keypair, encrypts mht_A and (-1)*tilde_mht_A, returns token for cloud.
    TokenA_Cloud Delegate();

    // Get the key pair (for passing pk to B and sk for decoding)
    const KeyPair& GetKeyPair() const { return kp_; }

    // Phase 4: Decode.
    // Decrypts 4 ciphertext vectors, recovers union via PRF + MHT indicator + BF queries.
    // Returns the union set, or empty vector on verification failure.
    std::vector<int64_t> Decode(const CloudResponse& response,
                                const TokenB_ToA& token_B);

    // Phase 5 (Algorithm 5): Insert (op='+') or delete (op='-') element x.
    // Returns true and fills out_req on success; false on duplicate insert
    // or non-existent delete. The caller is responsible for sending out_req
    // to the server. The element x must already lie in the PHF's domain.
    bool Update(char op, int64_t x, UpdateRequest& out_req);

    // Witness-refresh special case: delete the current witness and insert a
    // fresh one (incrementing the counter). Returns the two requests in
    // (delete, insert) order. The new witness must already be in the PHF
    // domain (the bench pre-includes a pool of reserve elements for this).
    bool RefreshWitness(int64_t new_witness, UpdateRequest& del_req,
                                              UpdateRequest& ins_req);

    int64_t GetWitness() const { return wit_A_val_; }

private:
    ProtocolConfig cfg_;
    FHEUtils& fhe_;
    const PHF& phf_;
    PRF prf_;
    MHT tilde_mht_;
    KeyPair kp_;

    NTL::ZZ wit_A_;
    int64_t wit_A_val_ = 0;

    // Bloom filters
    BloomFilter bf_A_;       // stores elements of X_A
    BloomFilter mask_bf_A_;  // stores masked encodings r*x+s
    BloomFilter xmask_bf_A_; // stores cross-masked encodings s*x+r

    // Plaintext modulus for FHE (from params.h)
    static constexpr int64_t P = static_cast<int64_t>(PLAINTEXT_MODULUS);

    static int64_t ModP(int64_t a);
    static int64_t MulModP(int64_t a, int64_t b);
    static int64_t ModInverse(int64_t a);
};

// ============================================================
// ClientB
// ============================================================
class ClientB {
public:
    ClientB(const ProtocolConfig& cfg, FHEUtils& fhe, const PHF& phf);

    // Compute the witness value wit_B = PRF(K3, ctr) mod P.
    NTL::ZZ ComputeWitness();

    // Phase 1: Encode input set X_B.
    // ComputeWitness() must have been called first.
    // Returns bar_e_B (the permuted masked encoding vector of length m).
    std::vector<int64_t> Encode(const std::vector<int64_t>& X_B);

    // Phase 2: Delegate.
    // Uses pk from A to encrypt mht_B, r_B, s_B.
    TokenB_Cloud Delegate(const KeyPair& pk_from_A);

    // Get the token for A (xmask_bf_B, wit_B)
    TokenB_ToA GetTokenForA() const;

    // Phase 5 (Algorithm 5): same as ClientA::Update but tagged as B.
    bool Update(char op, int64_t x, UpdateRequest& out_req);

    bool RefreshWitness(int64_t new_witness, UpdateRequest& del_req,
                                              UpdateRequest& ins_req);

    int64_t GetWitness() const { return wit_B_val_; }

private:
    ProtocolConfig cfg_;
    FHEUtils& fhe_;
    const PHF& phf_;
    PRF prf_;
    MHT tilde_mht_;

    NTL::ZZ wit_B_;
    int64_t wit_B_val_ = 0;

    // Bloom filters
    BloomFilter bf_B_;
    BloomFilter mask_bf_B_;
    BloomFilter xmask_bf_B_;

    static constexpr int64_t P = static_cast<int64_t>(PLAINTEXT_MODULUS);

    static int64_t ModP(int64_t a);
    static int64_t MulModP(int64_t a, int64_t b);
};

} // namespace vdpsu

#endif // VDPSU_CLIENT_H

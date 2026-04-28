#ifndef VDPSU_CLOUD_H
#define VDPSU_CLOUD_H

#include "params.h"
#include "fhe_utils.h"
#include "client.h"

#include <vector>

namespace vdpsu {

class Cloud {
public:
    explicit Cloud(const ProtocolConfig& cfg, FHEUtils& fhe);

    // Compute phase: given tokens from A and B, plus their masked encodings,
    // perform the homomorphic computation and return the response.
    //
    // token_A: (ct_mht_A, ct_neg_tilde_mht_A)
    // token_B: (ct_mht_B, ct_r_B, ct_s_B)
    // bar_e_A: plaintext masked encoding vector from A (length m)
    // bar_e_B: plaintext masked encoding vector from B (length m)
    //
    // Returns CloudResponse containing:
    //   ct_bar_y     = bar_e_A * [mht_A] + bar_e_B * [vecd]
    //   ct_vecd      = [mht_B] * [delta]
    //   ct_xmask_r_B = [r_B] * [vecd]
    //   ct_xmask_s_B = [s_B] * [vecd]
    // where delta = [mht_B] + [(-1)*tilde_mht_A]
    CloudResponse Compute(const TokenA_Cloud& token_A,
                          const TokenB_Cloud& token_B,
                          const std::vector<int64_t>& bar_e_A,
                          const std::vector<int64_t>& bar_e_B);

    // Phase 5 (Algorithm 5, server side): apply an in-place update to the
    // stored masked vector for the requesting party.
    void ProcessUpdate(const UpdateRequest& req,
                       std::vector<int64_t>& bar_e_party) const;

private:
    ProtocolConfig cfg_;
    FHEUtils& fhe_;
};

} // namespace vdpsu

#endif // VDPSU_CLOUD_H

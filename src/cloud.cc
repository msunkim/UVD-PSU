#include "cloud.h"

namespace vdpsu {

Cloud::Cloud(const ProtocolConfig& cfg, FHEUtils& fhe)
    : cfg_(cfg), fhe_(fhe) {}

CloudResponse Cloud::Compute(const TokenA_Cloud& token_A,
                              const TokenB_Cloud& token_B,
                              const std::vector<int64_t>& bar_e_A,
                              const std::vector<int64_t>& bar_e_B) {
    // delta = [mht_B] + [(-1)*tilde_mht_A]
    CTVector ct_delta = fhe_.HAdd(token_B.ct_mht_B, token_A.ct_neg_tilde_mht_A);

    // [vecd] = [mht_B] * [delta]  (depth-1 ct-ct multiplication)
    CTVector ct_vecd = fhe_.HMul(token_B.ct_mht_B, ct_delta);

    // [xmask_r_B] = [r_B] * [vecd]  (depth-1 ct-ct multiplication)
    CTVector ct_xmask_r_B = fhe_.HMul(token_B.ct_r_B, ct_vecd);

    // [xmask_s_B] = [s_B] * [vecd]  (depth-1 ct-ct multiplication)
    CTVector ct_xmask_s_B = fhe_.HMul(token_B.ct_s_B, ct_vecd);

    // [bar_y] = bar_e_A * [mht_A] + bar_e_B * [vecd]
    // Both are plaintext-ciphertext multiplications followed by addition
    CTVector ct_term_A = fhe_.PTMul(bar_e_A, token_A.ct_mht_A);
    CTVector ct_term_B = fhe_.PTMul(bar_e_B, ct_vecd);
    CTVector ct_bar_y = fhe_.HAdd(ct_term_A, ct_term_B);

    CloudResponse response;
    response.ct_bar_y = std::move(ct_bar_y);
    response.ct_vecd = std::move(ct_vecd);
    response.ct_xmask_r_B = std::move(ct_xmask_r_B);
    response.ct_xmask_s_B = std::move(ct_xmask_s_B);
    return response;
}

void Cloud::ProcessUpdate(const UpdateRequest& req,
                           std::vector<int64_t>& bar_e_party) const {
    if (req.pos < bar_e_party.size()) {
        bar_e_party[req.pos] = req.mask_x;
    }
}

} // namespace vdpsu

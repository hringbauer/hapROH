import numpy as np
from numba import njit

# Row meaning in transition probability matrix
STAY_OUT, ENTER_ROH, LEAVE_ROH, STAY_ROH, JUMP_ROH = (0,1,2,3,4)

# ---------------------------------------------------------------
# Numba version (from Claude)
# ---------------------------------------------------------------

@njit(cache=True, fastmath=True)
def _forward_backward_numba(ref_panel, proba_e, proba_t, nb_ref, nb_snp, nb_samples):
    post_pb = np.empty((nb_ref + 1, nb_snp, nb_samples))

    # ---------------------------------------------------------------
    # Init first SNP
    # ---------------------------------------------------------------
    for s in range(nb_samples):
        total = proba_e[2, 0, s]
        post_pb[0, 0, s] = total
        for r in range(nb_ref):
            v = proba_e[ref_panel[0, r], 0, s]
            post_pb[r + 1, 0, s] = v
            total += v
        inv_total = 1.0 / total
        for r in range(nb_ref + 1):
            post_pb[r, 0, s] *= inv_total

    # ---------------------------------------------------------------
    # Forward pass
    # ---------------------------------------------------------------
    for i in range(nb_snp - 1):
        pt_stay_out = proba_t[STAY_OUT, i]
        pt_enter = proba_t[ENTER_ROH, i]
        pt_leave = proba_t[LEAVE_ROH, i]
        pt_jump = proba_t[JUMP_ROH, i]
        coef = proba_t[STAY_ROH, i] - pt_jump

        for s in range(nb_samples):
            p0 = post_pb[0, i, s]
            sum_roh = 1.0 - p0

            val0 = (pt_stay_out * p0 + pt_leave * sum_roh) * proba_e[2, i, s]
            post_pb[0, i + 1, s] = val0
            total = val0

            base = pt_enter * p0 + pt_jump * sum_roh
            for r in range(nb_ref):
                emis = proba_e[ref_panel[i, r], i, s]
                v = (base + coef * post_pb[r + 1, i, s]) * emis
                post_pb[r + 1, i + 1, s] = v
                total += v

            inv_total = 1.0 / total
            for r in range(nb_ref + 1):
                post_pb[r, i + 1, s] *= inv_total

    # ---------------------------------------------------------------
    # Backward pass, combined in-place with the forward result
    # ---------------------------------------------------------------
    prev_bwd = np.full((nb_ref + 1, nb_samples), 1.0 / (nb_ref + 1))
    cur_bwd = np.empty((nb_ref + 1, nb_samples))
    terms = np.empty(nb_ref)  # scratch, allocated once, reused every (i, s)

    for i in range(nb_snp - 1, 0, -1):
        pt_stay_out = proba_t[STAY_OUT, i - 1]
        pt_enter = proba_t[ENTER_ROH, i - 1]
        pt_leave = proba_t[LEAVE_ROH, i - 1]
        pt_jump = proba_t[JUMP_ROH, i - 1]
        coef = proba_t[STAY_ROH, i - 1] - pt_jump

        for s in range(nb_samples):
            not_roh = proba_e[2, i, s] * prev_bwd[0, s]

            sum_roh = 0.0
            for r in range(nb_ref):
                emis = proba_e[ref_panel[i, r], i, s]
                t = emis * prev_bwd[r + 1, s]
                terms[r] = t
                sum_roh += t

            c0 = pt_stay_out * not_roh + pt_enter * sum_roh
            cur_bwd[0, s] = c0
            total = c0

            for r in range(nb_ref):
                v = coef * terms[r] + pt_leave * not_roh + pt_jump * sum_roh
                cur_bwd[r + 1, s] = v
                total += v

            inv_total = 1.0 / total
            for r in range(nb_ref + 1):
                cur_bwd[r, s] *= inv_total
                post_pb[r, i - 1, s] *= cur_bwd[r, s]

        prev_bwd, cur_bwd = cur_bwd, prev_bwd

    # ---------------------------------------------------------------
    # Final normalisation
    # ---------------------------------------------------------------
    for i in range(nb_snp):
        for s in range(nb_samples):
            total = 0.0
            for r in range(nb_ref + 1):
                total += post_pb[r, i, s]
            inv_total = 1.0 / total
            for r in range(nb_ref + 1):
                post_pb[r, i, s] *= inv_total

    return post_pb
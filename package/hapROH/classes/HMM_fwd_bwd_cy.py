"""
Pure-Python-style Cython implementation of the forward-backward recursion.

This file is valid, importable Python as-is (the `cython.*` annotations are
no-ops when not compiled). When run through `cythonize()`, the type
annotations let Cython generate a fast C implementation with no per-element
Python object overhead.
"""

import cython
import numpy as np

# Order of the states in TransitionProba (kept identical to forward_backward.py)
STAY_OUT: cython.int = 0
ENTER_ROH: cython.int = 1
LEAVE_ROH: cython.int = 2
STAY_ROH: cython.int = 3
JUMP_ROH: cython.int = 4


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _forward_backward_cython(
    ref_panel: cython.uchar[:, :],
    proba_e: cython.double[:, :, :],
    proba_t: cython.double[:, :],
    nb_ref: cython.int,
    nb_snp: cython.int,
    nb_samples: cython.int,
):
    """
    Same as `_forward_backward` / `_forward_backward_numba`, but compiled
    via Cython. Args:
        ref_panel: (nb_snp, nb_ref) uint8 memoryview, 0 (REF) or 1 (ALT)
        proba_e:   (3, nb_snp, nb_samples) float64 memoryview
        proba_t:   (5, nb_snp-1) float64 memoryview
    Returns:
        post_pb: np.ndarray of shape (nb_ref+1, nb_snp, nb_samples)
    """
    post_pb_arr = np.empty((nb_ref + 1, nb_snp, nb_samples), dtype=np.float64)
    post_pb: cython.double[:, :, :] = post_pb_arr

    # --- typed locals: these compile to C variables, not Python objects ---
    i: cython.int
    r: cython.int
    s: cython.int
    total: cython.double
    inv_total: cython.double
    p0: cython.double
    sum_roh: cython.double
    val0: cython.double
    base: cython.double
    v: cython.double
    emis: cython.double
    coef: cython.double
    pt_stay_out: cython.double
    pt_enter: cython.double
    pt_leave: cython.double
    pt_jump: cython.double
    not_roh: cython.double
    c0: cython.double
    t: cython.double

    # ---- Init first SNP ----
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

    # ---- Forward pass ----
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

    # ---- Backward pass, combined in-place with forward result ----
    prev_bwd_arr = np.full((nb_ref + 1, nb_samples), 1.0 / (nb_ref + 1), dtype=np.float64)
    cur_bwd_arr = np.empty((nb_ref + 1, nb_samples), dtype=np.float64)
    terms_arr = np.empty(nb_ref, dtype=np.float64)

    prev_bwd: cython.double[:, :] = prev_bwd_arr
    cur_bwd: cython.double[:, :] = cur_bwd_arr
    terms: cython.double[:] = terms_arr

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

        # swap buffers (plain Python-level tuple assignment, cheap either way)
        prev_bwd, cur_bwd = cur_bwd, prev_bwd

    # ---- Final normalisation ----
    for i in range(nb_snp):
        for s in range(nb_samples):
            total = 0.0
            for r in range(nb_ref + 1):
                total += post_pb[r, i, s]
            inv_total = 1.0 / total
            for r in range(nb_ref + 1):
                post_pb[r, i, s] *= inv_total

    return post_pb_arr
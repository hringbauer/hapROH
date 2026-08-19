import logging

import numpy as np
from numba import njit

from .genomicData import GenomicData, DataType
from .transitionProba import TransitionProba, get_transi_proba
from .emissionProba import EmissionProba, get_emi_proba
from .HMM_fwd_bwd_cy import _forward_backward_cython

from hapROH.utils.miscellanious import print_memory_usage

logger = logging.getLogger(__name__)

################################
# Forward backward algorithms
################################
"""
Note: data used
input:
    ref_panel: (nb_snp, nb_ref), uint8 containing 0 (REF) or 1 (ALT)
    proba_e: (3, nb_snp, nb_samples), float, with 3 states: ROH_REF, ROH_ALT, no_ROH
    proba_t: (5, nb_snp-1), float, with 5 states (see below)

output:
    post_pb: (nb_ref+1, nb_snp, nb_samples), float, with n_ref+1 states:
        no_ROH in row 0 and ROH_with_ref_i in row i+1
"""
# Order of the states in TransitionProba
STAY_OUT, ENTER_ROH, LEAVE_ROH, STAY_ROH, JUMP_ROH = (0,1,2,3,4)

def _forward_backward(ref_panel:np.ndarray, proba_e:EmissionProba, proba_t:TransitionProba, nb_ref:int, nb_snp:int, nb_samples:int) -> np.ndarray:
    """
    Compute the posterior probability of each state at each locus,
    using the standard forward-backward algorithm.
    Values are normalised at each step to avoid numerical issues (because probas converge to 0)
    Returns:
        post_pb: np.ndarray of shape (nb_ref+1, nb_snp, nb_samples)
    """
    nb_snp, nb_ref = ref_panel.shape
    _,_, nb_samples = proba_e.shape
    post_pb = np.empty((nb_ref+1, nb_snp, nb_samples), dtype=float)

    ### Initialize first SNP
    post_pb[0, 0] = proba_e[2, 0]                 # not ROH
    post_pb[1:, 0] = proba_e[ref_panel[0], 0]     # ROH with one of the n_ref possibilities
    post_pb[:, 0] /= post_pb[:, 0].sum(axis=0)    # rescale to 1

    ### Preallocated scratch buffers, reused every iteration
    emission_roh = np.empty((nb_ref, nb_samples), dtype=float)
    term = np.empty((nb_ref, nb_samples), dtype=float)

    ### Forward algorithm
    logger.debug("Starting forward computation")
    for i in range(0, nb_snp-1):
        sum_roh = 1 - post_pb[0, i]         # = post_pb[1:].sum(axis=0) because of normalisation

        ### Non ROH state
        post_pb[0, i+1] = proba_t[STAY_OUT, i] * post_pb[0, i]
        post_pb[0, i+1] += proba_t[LEAVE_ROH, i] * sum_roh
        post_pb[0, i+1] *= proba_e[2, i]

        ### ROH states
        # Emission probas gathered into preallocated buffer
        np.take(proba_e[:, i], ref_panel[i], axis=0, out=emission_roh)

        # (STAY_ROH - JUMP_ROH) * post_pb[1:, i]  -> reuse preallocated buffer `term`
        coef = proba_t[STAY_ROH, i] - proba_t[JUMP_ROH, i]
        np.multiply(post_pb[1:, i], coef, out=term)

        post_pb[1:, i+1] = proba_t[ENTER_ROH, i] * post_pb[0, i]
        post_pb[1:, i+1] += proba_t[JUMP_ROH, i] * sum_roh
        post_pb[1:, i+1] += term
        post_pb[1:, i+1] *= emission_roh

        post_pb[:, i+1] /= post_pb[:, i+1].sum(axis=0)
    logger.debug("Done forward computation")
    print_memory_usage(logger)

    ### Backward algorithm
    # note: bwd is directly combined with fwd. only prev_bwd is stored, as needed for computation
    logger.debug("Starting backward computation")
    prev_bwd = np.ones((nb_ref+1, nb_samples), dtype=float) / (nb_ref+1)
    cur_bwd =  np.empty((nb_ref+1, nb_samples), dtype=float)
    for i in range(nb_snp-1, 0, -1):
        # Gather emission probas into preallocated buffer
        np.take(proba_e[:, i], ref_panel[i], axis=0, out=emission_roh)

        # emission_roh * prev_bwd[1:]  -> preallocated buffer `term`
        np.multiply(emission_roh, prev_bwd[1:], out=term)
        sum_roh = term.sum(axis=0)
        not_roh = proba_e[2, i] * prev_bwd[0]

        ### Non ROH state
        cur_bwd[0] = proba_t[STAY_OUT, i-1] * not_roh
        cur_bwd[0] += proba_t[ENTER_ROH, i-1] * sum_roh

        ### ROH states
        coef = proba_t[STAY_ROH, i - 1] - proba_t[JUMP_ROH, i - 1]
        np.multiply(term, coef, out=cur_bwd[1:])
        cur_bwd[1:] += proba_t[LEAVE_ROH, i-1] * not_roh
        cur_bwd[1:] += proba_t[JUMP_ROH, i-1] * sum_roh

        cur_bwd /= cur_bwd.sum(axis=0)
        post_pb[:, i-1] *= cur_bwd
        prev_bwd, cur_bwd = cur_bwd, prev_bwd
    logger.debug("Done backward computation")
    print_memory_usage(logger)

    post_pb /= post_pb.sum(axis=0, keepdims=True)

    return post_pb

@njit(cache=True, fastmath=True)
def _forward_backward_numba(ref_panel, proba_e, proba_t, nb_ref, nb_snp, nb_samples):
    """Same as `_forward_backward` but with numba implementation."""
    post_pb = np.empty((nb_ref + 1, nb_snp, nb_samples))

    # Init first SNP
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

    # Forward pass
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

    # Backward pass, combined in-place with the forward result
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

    # Final normalisation
    for i in range(nb_snp):
        for s in range(nb_samples):
            total = 0.0
            for r in range(nb_ref + 1):
                total += post_pb[r, i, s]
            inv_total = 1.0 / total
            for r in range(nb_ref + 1):
                post_pb[r, i, s] *= inv_total

    return post_pb

################################
# HMM class
################################

class HMM():
    ### Data
    ref_panel:np.ndarray            # shape (nb_snp, n_ref), dtype 0|1

    ### Transition parameters
    r_in: float             # Rate of jump into ROH state (per Morgan)
    r_out: float            # Rate of jump out of ROH state (per Morgan)
    r_jump: float           # Rate of jump between two ROH states (per Morgan)

    ### Emission parameters
    error_rate: float       # Sequencing error rate.

    ### Probabilities
    proba_t: TransitionProba    # shape (5, nb_snp-1, nb_samples), dtype float, contains following proba: STAY_OUT, ENTER_ROH, LEAVE_ROH, STAY_ROH, JUMP_ROH
    proba_e: EmissionProba      # shape (3, nb_snp, nb_samples), dtype float, describe proba between the following states: ROH_REF, ROH_ALT, no_ROH

    def __init__(self, sample_data: GenomicData, ref_data:GenomicData, r_map:np.ndarray,
                    r_in:float, r_out:float, r_jump: float,
                    error_rate:float
                ) -> None:
        """Initialize HMM by computing emission and transition probabilities"""
        if ref_data.datatype != DataType.GT:
            raise ValueError(f"Expecte reference pannel to contain GT, but contains datatype {ref_data.datatype}")
        self.ref_panel = ref_data.data.reshape(ref_data.data.shape[0], -1)   # (nb_snp, nb_samples, 2) -> (nb_snp, 2*nb_samples)
    
        self.r_in = r_in
        self.r_out = r_out
        self.r_jump = r_jump

        self.error_rate = error_rate

        print_memory_usage(logger)
        logger.debug("Computing transition probabilities")
        self.proba_t = get_transi_proba(r_in, r_out, r_jump, self.ref_panel.shape[1], r_map)     # shape (5, nb_snp-1)
        print_memory_usage(logger)

        logger.debug("Computing emission probabilities")
        allele_freq = self.ref_panel.mean(axis=1)
        self.proba_e = get_emi_proba(sample_data, allele_freq, error_rate)                  # shape (3, nb_snp, nb_samples)
        print_memory_usage(logger)

    def calc_posterior_proba(self, backend:Literal["python", "numba", "cython"]="cython") -> np.ndarray:
        """
        Compute the posterior probability of each state at each locus.
        Args:
            backend: Which backend to the for the main computation algorithm.
        Returns:
            post_pb: np.ndarray of shape (nb_ref+1, nb_snp, nb_samples)
        """
        logger.debug("Computing posterior probabilities")
        nb_snp, nb_ref = self.ref_panel.shape
        _, _, nb_samples = self.proba_e.shape
        ref_panel = np.ascontiguousarray(self.ref_panel, dtype=np.uint8)
        proba_e = np.ascontiguousarray(self.proba_e, dtype=np.float64)
        proba_t = np.ascontiguousarray(self.proba_t, dtype=np.float64)
        match backend:
            case "numba":
                return _forward_backward_numba(ref_panel, proba_e, proba_t, nb_ref, nb_snp, nb_samples)
            case "cython":
                return _forward_backward_cython(ref_panel, proba_e, proba_t, nb_ref, nb_snp, nb_samples)
            case "python":
                return _forward_backward(ref_panel, proba_e, proba_t, nb_ref, nb_snp, nb_samples)
        raise ValueError(f"Unknown backend ({backend}). Available implementations are: 'python', 'numba' or 'cython'")

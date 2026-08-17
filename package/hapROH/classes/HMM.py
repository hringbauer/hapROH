import logging

import numpy as np

from .genomicData import GenomicData, DataType
from .transitionProba import TransitionProba, get_transi_proba
from .emissionProba import EmissionProba, get_emi_proba
from .fwd_bwd import _forward_backward_numba

from hapROH.utils.miscellanious import print_memory_usage

logger = logging.getLogger(__name__)

# Row meaning in transition probability
STAY_OUT, ENTER_ROH, LEAVE_ROH, STAY_ROH, JUMP_ROH = (0,1,2,3,4)

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

    def calc_posterior_proba(self) -> np.ndarray:
        """
        Compute the posterior probability of each state at each locus.
        Values are normalised at each step to avoid numerical issues (because probas converge to 0)
        Returns:
            post_pb: np.ndarray of shape (nb_ref+1, nb_snp, nb_samples)
        """
        logger.info("Computing posterior probabilities")
        nb_snp, nb_ref = self.ref_panel.shape
        _,_, nb_samples = self.proba_e.shape
        post_pb = np.empty((nb_ref+1, nb_snp, nb_samples), dtype=float)

        ### Initialize first SNP
        post_pb[0, 0] = self.proba_e[2, 0]                      # not ROH
        post_pb[1:, 0] = self.proba_e[self.ref_panel[0], 0]     # ROH with one of the n_ref possibilities
        post_pb[:, 0] /= post_pb[:, 0].sum(axis=0)              # rescale to 1

        ### Preallocated scratch buffers, reused every iteration
        emission_roh = np.empty((nb_ref, nb_samples), dtype=float)
        term = np.empty((nb_ref, nb_samples), dtype=float)

        ### Forward algorithm
        logger.debug("Starting forward computation")
        for i in range(0, nb_snp-1):
            sum_roh = 1 - post_pb[0, i]         # = post_pb[1:].sum(axis=0) because of normalisation

            ### Non ROH state
            post_pb[0, i+1] = self.proba_t[STAY_OUT, i] * post_pb[0, i]
            post_pb[0, i+1] += self.proba_t[LEAVE_ROH, i] * sum_roh
            post_pb[0, i+1] *= self.proba_e[2, i]

            ### ROH states
            # Emission probas gathered into preallocated buffer
            np.take(self.proba_e[:, i], self.ref_panel[i], axis=0, out=emission_roh)

            # (STAY_ROH - JUMP_ROH) * post_pb[1:, i]  -> reuse preallocated buffer `term`
            coef = self.proba_t[STAY_ROH, i] - self.proba_t[JUMP_ROH, i]
            np.multiply(post_pb[1:, i], coef, out=term)

            post_pb[1:, i+1] = self.proba_t[ENTER_ROH, i] * post_pb[0, i]
            post_pb[1:, i+1] += self.proba_t[JUMP_ROH, i] * sum_roh
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
            np.take(self.proba_e[:, i], self.ref_panel[i], axis=0, out=emission_roh)

            # emission_roh * prev_bwd[1:]  -> preallocated buffer `term`
            np.multiply(emission_roh, prev_bwd[1:], out=term)
            sum_roh = term.sum(axis=0)
            not_roh = self.proba_e[2, i] * prev_bwd[0]

            ### Non ROH state
            cur_bwd[0] = self.proba_t[STAY_OUT, i-1] * not_roh
            cur_bwd[0] += self.proba_t[ENTER_ROH, i-1] * sum_roh

            ### ROH states
            coef = self.proba_t[STAY_ROH, i - 1] - self.proba_t[JUMP_ROH, i - 1]
            np.multiply(term, coef, out=cur_bwd[1:])
            cur_bwd[1:] += self.proba_t[LEAVE_ROH, i-1] * not_roh
            cur_bwd[1:] += self.proba_t[JUMP_ROH, i-1] * sum_roh

            cur_bwd /= cur_bwd.sum(axis=0)
            post_pb[:, i-1] *= cur_bwd
            prev_bwd, cur_bwd = cur_bwd, prev_bwd
        logger.debug("Done backward computation")
        print_memory_usage(logger)

        post_pb /= post_pb.sum(axis=0, keepdims=True)

        return post_pb

    def calc_posterior_proba_numba(self) -> np.ndarray:
        """
        Compute the posterior probability of each state at each locus.
        Values are normalised at each step to avoid numerical issues (because probas converge to 0).
        Same as `calc_posterior_proba` but with numba compilation.
        Returns:
            post_pb: np.ndarray of shape (nb_ref+1, nb_snp, nb_samples)
        """
        nb_snp, nb_ref = self.ref_panel.shape
        _, _, nb_samples = self.proba_e.shape
    
        # ref_panel must be a plain int array for numba indexing
        ref_panel = np.ascontiguousarray(self.ref_panel)
        proba_e = np.ascontiguousarray(self.proba_e, dtype=np.float64)
        proba_t = np.ascontiguousarray(self.proba_t, dtype=np.float64)
    
        return _forward_backward_numba(ref_panel, proba_e, proba_t, nb_ref, nb_snp, nb_samples)

import logging

import numpy as np

from .genomicData import GenomicData
from .transitionProba import TransitionProba, get_transi_proba
from .emissionProba import EmissionProba, get_emi_proba

from hapROH.utils.miscellanious import print_memory_usage

logger = logging.getLogger(__name__)

def test():
    pass

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
    proba_t: TransitionProba    # shape (5, nb_snp-1, nb_samples), dtype float, contains following proba: stay_out, enter_ROH, leave_ROH, stay_ROH, jump_ROH
    proba_e: EmissionProba      # shape (3, nb_snp, nb_samples), dtype float, describe proba between the following states: ROH_REF, ROH_ALT, no_ROH

    def __init__(self, sample_data: GenomicData, ref_panel:np.ndarray, r_map:np.ndarray,
                    r_in:float, r_out:float, r_jump: float,
                    error_rate:float
                ) -> None:
        """Initialize HMM by computing emission and transition probabilities"""
        self.ref_panel = ref_panel                                                          # shape (nb_snp, nb_ref)

        self.r_in = r_in
        self.r_out = r_out
        self.r_jump = r_jump

        self.error_rate = error_rate

        print_memory_usage(logger)
        logger.debug("Computing transition probabilities")
        self.proba_t = get_transi_proba(r_in, r_out, r_jump, ref_panel.shape[1], r_map)     # shape (5, nb_snp-1)
        print_memory_usage(logger)

        logger.debug("Computing emission probabilities")
        allele_freq = ref_panel.mean(axis=1)
        self.proba_e = get_emi_proba(sample_data, allele_freq, error_rate)                  # shape (3, nb_snp, nb_samples)
        print_memory_usage(logger)

    def calc_posterior_proba(self) -> np.ndarray:
        """
        Compute the posterior probability of each state at each locus.
        Values are normalised at each step to avoid numerical issues (proba converging to 0)
        Returns:
            post_pb: np.ndarray of shape (nb_ref+1, nb_snp, nb_samples)
        """
        logger.info("Computing posterior probabilities")
        nb_snp, nb_ref = self.ref_panel.shape
        _,_, nb_samples = self.proba_e.shape
        post_pb = np.empty((nb_ref+1, nb_snp, nb_samples), dtype=float)

        # Row order in self.proba_t
        stay_out, enter_ROH, leave_ROH, stay_ROH, jump_ROH = (0,1,2,3,4)

        ### Initialize first SNP
        post_pb[0, 0] = self.proba_e[2, 0]                      # not ROH
        post_pb[1:, 0] = self.proba_e[self.ref_panel[0], 0]     # ROH with one of the n_ref possibilities
        post_pb[:, 0] /= post_pb[:, 0].sum(axis=0)              # rescale to 1

        ### Preallocated scratch buffers, reused every iteration
        emission_roh = np.empty((nb_ref, nb_samples), dtype=float)
        term = np.empty((nb_ref, nb_samples), dtype=float)

        ### Forward algorithm
        logger.debug(f"(nb_ref, nb_snp, nb_samples: {(nb_ref, nb_snp, nb_samples)}")
        print_memory_usage(logger)
        logger.debug("Starting forward computation")
        for i in range(0, nb_snp-1):
            sum_roh = 1 - post_pb[0, i]         # = post_pb[1:].sum(axis=0) because of normalisation

            ### Non ROH state
            post_pb[0, i+1] = self.proba_t[stay_out, i] * post_pb[0, i]
            post_pb[0, i+1] += self.proba_t[leave_ROH, i] * sum_roh
            post_pb[0, i+1] *= self.proba_e[2, i]

            ### ROH states
            # Emission probas gathered into preallocated buffer
            np.take(self.proba_e[:, i], self.ref_panel[i], axis=0, out=emission_roh)

            # (stay_ROH - jump_ROH) * post_pb[1:, i]  -> reuse preallocated buffer `term`
            coef = self.proba_t[stay_ROH, i] - self.proba_t[jump_ROH, i]
            np.multiply(post_pb[1:, i], coef, out=term)

            post_pb[1:, i+1] = self.proba_t[enter_ROH, i] * post_pb[0, i]
            post_pb[1:, i+1] += self.proba_t[jump_ROH, i] * sum_roh
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
            cur_bwd[0] = self.proba_t[stay_out, i-1] * not_roh
            cur_bwd[0] += self.proba_t[enter_ROH, i-1] * sum_roh

            ### ROH states
            coef = self.proba_t[stay_ROH, i - 1] - self.proba_t[jump_ROH, i - 1]
            np.multiply(term, coef, out=cur_bwd[1:])
            cur_bwd[1:] += self.proba_t[leave_ROH, i-1] * not_roh
            cur_bwd[1:] += self.proba_t[jump_ROH, i-1] * sum_roh

            cur_bwd /= cur_bwd.sum(axis=0)
            post_pb[:, i-1] *= cur_bwd
            prev_bwd, cur_bwd = cur_bwd, prev_bwd
        logger.debug("Done backward computation")
        print_memory_usage(logger)

        post_pb /= post_pb.sum(axis=0, keepdims=True)

        return post_pb
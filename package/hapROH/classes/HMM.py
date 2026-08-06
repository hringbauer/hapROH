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
    proba_e: EmissionProba      # shape (3, nb_snp, nb_samples), dtype float, describe proba between the following states: no_ROH, ROH_REF, ROH_ALT

    def __init__(self, sample_data: GenomicData, ref_panel:np.ndarray, r_map:np.ndarray,
                    r_in:float, r_out:float, r_jump: float,
                    error_rate:float
                ) -> None:
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
        Compute the posterior probability of each state at each locus
        Returns:
            post_pb: np.ndarray of shape (nb_ref+1, nb_snp, nb_samples)
        """
        logger.info("Computing posterior probabilities")
        nb_snp, nb_ref = self.ref_panel.shape
        _,_, nb_samples = self.proba_e.shape
        fwd = np.empty((nb_ref+1, nb_snp, nb_samples), dtype=float)
        bwd = np.empty((nb_ref+1, nb_snp, nb_samples), dtype=float)

        # Row order in self.proba_t
        stay_out, enter_ROH, leave_ROH, stay_ROH, jump_ROH = (0,1,2,3,4)

        ### Initialize first/last SNP
        fwd[0, 0] = self.proba_e[0, 0]                        # not ROH
        fwd[1:, 0] = self.proba_e[self.ref_panel[0, :], 0]    # ROH with one of the n_ref possibilities
        bwd[:, -1] = 1

        ### Forward algorithm
        logger.debug("Starting forward computation")
        for i in range(0, nb_snp-1):
            fwd[0, i+1] = self.proba_e[0, i] * (self.proba_t[stay_out, i] * fwd[0, i]
                                                 + self.proba_t[leave_ROH, i] * fwd[1:, i].sum(axis=0))
            fwd[1:, i+1] = self.proba_e[self.ref_panel[i, :], i] * (
                  self.proba_t[enter_ROH, i] * fwd[0, i]
                + self.proba_t[jump_ROH, i] * fwd[1:, i].sum(axis=0)
                + (self.proba_t[stay_ROH, i] - self.proba_t[jump_ROH, i] ) * fwd[1:, i]
            )
        logger.debug("Done forward computation")

        ### Backward algorithm
        logger.debug("Starting backward computation")
        for i in range(nb_snp-1, 0, -1):
            bwd[0, i-1] = self.proba_t[stay_out, i-1] * (self.proba_e[0, i] * bwd[0, i]) \
                    + self.proba_t[enter_ROH, i-1] * (self.proba_e[self.ref_panel[i, :], i] * bwd[1:, i]).sum(axis=0)
            bwd[1:, i-1] = self.proba_t[leave_ROH, i-1] * (self.proba_e[0, i] * bwd[0, i]) \
                    + self.proba_t[jump_ROH, i-1] * (self.proba_e[self.ref_panel[i, :], i] * bwd[1:, i]).sum(axis=0) \
                    + (self.proba_t[stay_ROH, i] - self.proba_t[jump_ROH, i] ) * (self.proba_e[1:, i] * bwd[1:, i])
        logger.debug("Done backward computation")

        # TODO/note: in Harald's/Yilei's implementation, fwd and bwd are scaled/taken into log
        # -> necessary for numeric precision because the value gets very small ?

        ### Combine 
        post_pb = fwd * bwd / fwd[:, -1].sum(axis=0)    # shape (nb_ref+1, nb_snp, nb_samples)

        return post_pb


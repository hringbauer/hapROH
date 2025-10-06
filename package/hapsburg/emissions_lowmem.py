"""
Class for calculating Emission Probabilities in a memory-efficient way.
Specifically, it computes emission probability for the non-ROH state and the ROH state (only two states, 0 or 1).
The original implementation computes the emission probability for all ROH states, 
which scales linearly with the number of haplotypes in the reference panel.
Contains Sub-Classes, as well as a factory Method.
@ Author: Yilei Huang, 2024, All rights reserved
"""

from hapsburg.emissions import Emissions
import numpy as np
from scipy.stats import binom


class Emissions(object):
    """Class for emission probabilities
    Has methods to return emission probabilities"""
    p = []  # Vector of alle frequencies [n_loci]
    ref_haps = []  # Reference Haplotypes [n_ref, n_loci/blocksize + overhang]
    e_rate = 1e-2  # Error rate of sequencing
    e_rate_ref = 1e-3  # The Probability of an error of copying from the reference
    overhang = 0  # Overhang when encoding genotypes in binary (blocks of 8 bits)
    blocksize = 8  # Blocksize for encoding genotypes in binary

    def __init__(self, ref_haps=[], blocksize=8, overhang=0):
        """Initialize Class"""
        if len(ref_haps) > 0:
            self.ref_haps = ref_haps
        self.blocksize = blocksize
        self.overhang = overhang

        # Calculate the allele frequencies
        self.p = self.calc_allelefreq_blockwise()
    

    def allele_freq_per_block(self, gts_one_block):
        """Calculate allele frequency per block.
        Return 1D array of length blocksize"""
        freq = np.empty(self.blocksize, dtype=float)
        for i in range(self.blocksize):
            freq[i] = np.mean(gts_one_block >> (self.blocksize -1 - i) & 1)
        return freq

    def calc_allelefreq_blockwise(self):
        """
        Calculate allele frequency for all markers in ref_haps, blockwise
        """
        allelefreq_blockwise = np.concatenate(np.apply_along_axis(self.allele_freq_per_block, 0, self.ref_haps).T, axis=0)
        if self.overhang > 0:
            allelefreq_blockwise = allelefreq_blockwise[:-(self.blocksize - self.overhang)]
        return allelefreq_blockwise

    def give_emission(self, ob_stat):
        """Gives the emission matrix of path of states
        Return emission matrix [3,l]"""
        raise NotImplementedError("Implement This in specific subclass.")

    def set_params(self, **kwargs):
        """Set the Parameters.
        Takes keyworded arguments"""
        for key, value in kwargs.items():
            setattr(self, key, value)


class Model_Emissions(Emissions):
    """Implements the haploid model Emission probabilities"""

    def give_emission(self, ob_stat):
        """Return the full emission Probability directly in Log Space.
        ob_stat: Observed Genotypes [2,l] (only use 1st row)"""

        n_loci = self.blocksize*np.shape(self.ref_haps)[1]
        if self.overhang > 0:
            n_loci -= (self.blocksize - self.overhang)
        ob_stat = ob_stat[0, :]  # Do ONLY first observed Haplotype
        assert(len(ob_stat) == n_loci)  # Sanity Check
        e_mat = np.zeros((3, n_loci), dtype=float)

        # Do the HW
        e_mat[0] = (ob_stat == 1) * (self.p * (1 - self.e_rate) + (1 - self.p) * self.e_rate) + \
                      (ob_stat == 0) * ((1 - self.p) * (1 - self.e_rate) + self.p * self.e_rate)
        # e_mat[1,:] emission probability when copied from reference allele
        # e_mat[2,:] emission probability when copied from alternative allele
        e_mat[1] = (ob_stat == 0) * (1 - self.e_rate) + (ob_stat == 1) * self.e_rate
        e_mat[2] = 1 - e_mat[1]
        assert(np.min(e_mat) >= 0)  # Sanity Check (In Normal Space Pr. >=0)
        assert(np.max(e_mat) <= 1)  # Sanity Check (In Normal Space Pr. <=1)
        return e_mat

class RC_Model_Emissions(Emissions):
    """Implements the Read Count model Emission probabilities."""

    def give_emission(self, ob_stat):
        assert(len(ob_stat) == 2)  # Sanity Check
        n_loci = self.blocksize*np.shape(self.ref_haps)[1]
        if self.overhang > 0:
            n_loci -= (self.blocksize - self.overhang)
        assert(len(ob_stat[0]) == n_loci)  # Sanity Check
        e_mat = np.zeros((3, n_loci), dtype=float)
        # Do the HW
        hw_prob = np.empty((n_loci, 3), dtype=float)
        hw_prob[:,0] = (1-self.p)**2
        hw_prob[:,1] = 2*self.p*(1-self.p)
        hw_prob[:,2] = self.p**2
        
        derived_read_prob = np.empty((1, 3), dtype=float)
        derived_read_prob[0, 0] = self.e_rate
        derived_read_prob[0, 1] = 0.5
        derived_read_prob[0, 2] = 1 - self.e_rate
        binom_pmf = binom.pmf(ob_stat[1].reshape(n_loci, 1), (ob_stat[0]+ob_stat[1]).reshape(n_loci,1), derived_read_prob) # shape (nloci, 3)
        e_mat[0] = np.sum(hw_prob*binom_pmf, axis=1)
        # Do the ROH state
        # e_mat[1,:] emission probability when copied from reference allele
        # e_mat[2,:] emission probability when copied from alternative allele
        genotype_prob = np.empty((2, 3), dtype=float)
        genotype_prob[0,0] = 1 - self.e_rate_ref
        genotype_prob[0,1] = self.e_rate_ref/2
        genotype_prob[0,2] = self.e_rate_ref/2
        genotype_prob[1,0] = self.e_rate_ref/2
        genotype_prob[1,1] = self.e_rate_ref/2
        genotype_prob[1,2] = 1 - self.e_rate_ref
        e_mat[1:] = genotype_prob@(binom_pmf.T)
        return e_mat

class Diploid_GT_Emissions(Emissions):
    """Implements the Diploid Genotype model Emission probabilities."""
    def give_emission(self, ob_stat):
        derived_gt = np.sum(ob_stat, axis=0)  # The Nr of of derived GT (0/1/2)
        e_rate = self.e_rate  # Error Rate (Here FLIP error)

        n_loci = self.blocksize*np.shape(self.ref_haps)[1]
        if self.overhang > 0:
            n_loci -= (self.blocksize - self.overhang)
        assert(len(ob_stat[0]) == n_loci)  # Sanity Check
        assert(len(derived_gt) == n_loci)  # Sanity Check

        e_mat = np.zeros((3, n_loci), dtype=float)

        ### HW State
        e_mat[0, :] = (1 - self.p) * (1 - self.p) * (derived_gt == 0) + \
            (1 - self.p) * self.p * 2 * (derived_gt == 1) + \
            self.p * self.p * (derived_gt == 2)
        ### Copying States
        e_mat[1:, :] = derived_gt == 0
        e_mat[2:, :] = derived_gt == 2
        e_mat = e_mat * (1 - self.e_rate) + (1 - e_mat) * self.e_rate / 2.0

        return e_mat

class RC_Model_Emissions_withContamination(Emissions):
    c = 0.0 # contamination rate
    pCon = [] # allele frequency of the contamination source

    def __init__(self, ref_haps, c, pCon, blocksize=8, overhang=0):
        super().__init__(ref_haps, blocksize, overhang)
        self.c = c
        self.pCon = pCon
    
    def give_emission(self, ob_stat):
        assert(len(ob_stat) == 2)
        n_loci = self.blocksize*np.shape(self.ref_haps)[1]
        if self.overhang > 0:
            n_loci -= (self.blocksize - self.overhang)
        assert(len(ob_stat[0]) == n_loci)  # Sanity Check
        e_mat = np.zeros((3, n_loci), dtype=float)

        c = self.c
        e_rate = self.e_rate
        pCon = self.pCon
        p_read = np.empty((n_loci, 3), dtype=float) # probability of sampling a derived read if the underlying genotype is 00,01,11, for each locus
        p_read[:,0] = (1-c)*e_rate + c*pCon*(1-e_rate) + c*(1-pCon)*e_rate
        p_read[:,1] = 0.5*(1-c) + c*pCon*(1-e_rate)
        p_read[:,2] = (1-c)*(1-e_rate) + c*pCon*(1-e_rate) + c*(1-pCon)*e_rate
        binom_pmf = binom.pmf(ob_stat[1].reshape(n_loci, 1), (ob_stat[0]+ob_stat[1]).reshape(n_loci,1), p_read) # shape: (n_loci, 3)

        hw_prob = np.empty((n_loci, 3), dtype=float)
        hw_prob[:,0] = (1-self.p)**2
        hw_prob[:,1] = 2*self.p*(1-self.p)
        hw_prob[:,2] = self.p**2
        e_mat = np.empty((3, n_loci), dtype=float)
        e_mat[0] = np.sum(hw_prob*binom_pmf, axis=1)
        # Do the ROH state
        genotype_prob = np.empty((2, 3), dtype=float)
        genotype_prob[0,0] = 1 - self.e_rate_ref
        genotype_prob[0,1] = 0.0
        genotype_prob[0,2] = self.e_rate_ref
        genotype_prob[1,0] = self.e_rate_ref
        genotype_prob[1,1] = 0.0
        genotype_prob[1,2] = 1 - self.e_rate_ref
        e_mat[1:] = genotype_prob@(binom_pmf.T)
        return e_mat

def load_emission_model_lowmem(ref_states, e_model="haploid", c=0.0, pCon=[], overhang=0):
    """Load the Emission Model"""
    if e_model == "haploid":
        e_obj = Model_Emissions(ref_states, overhang=overhang)
    elif e_model == "readcount":
        e_obj = RC_Model_Emissions(ref_states, overhang=overhang)
    elif e_model == "diploid_gt":
        e_obj = Diploid_GT_Emissions(ref_states, overhang=overhang)
    elif e_model == "readcount_contam":
        e_obj = RC_Model_Emissions_withContamination(ref_states, c, pCon, overhang=overhang)
    else:
        raise NotImplementedError("Emission Model not found!")
    return e_obj
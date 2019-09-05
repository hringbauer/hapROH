"""
Class for calculating Emission Probabilities.
Contains Sub-Classes, as well as factory Method.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import numpy as np
from scipy.stats import binom  # Binomial Likelihood

###############################
###############################


class Emissions(object):
    """Class for emission probabilities
    Has methods to return emission probabilities"""

    def give_emission_matrix(self, remember=False):
        """Return Emission Matrix"""
        raise NotImplementedError("Implement This in specific subclass.")

    def give_emission_state(self, ob_stat):
        """Gives the emission matrix of path of states"""
        raise NotImplementedError("Implement This in specific subclass.")

    def set_params(self, **kwargs):
        """Set the Parameters.
        Takes keyworded arguments"""
        for key, value in kwargs.items():
            setattr(self, key, value)


class Model_Emissions(Emissions):
    """Implements the haploid model Emission probabilities"""
    p = []  # Vector of alle frequencies [n_loci]
    ref_haps = []  # Reference Haplotypes [n_ref, n_loci]
    e_mat = []  # Full Emission Matrix [n_ref+1, n_loci, 2]
    e_rate = 1e-3  # The Probability of an error (e-10: Default)

    def __init__(self, ref_haps=[]):
        """Initialize Class"""
        if len(ref_haps) > 0:
            self.ref_haps = ref_haps

        # Calculate the allele frequencies
        self.p = np.mean(self.ref_haps, axis=0)

    def give_emission_matrix(self, remember=False):
        """Return full Emission Matrix.
        dtype: Precision of the returned Matrix"""
        n_loci = np.shape(self.ref_haps)[1]
        n_ref = np.shape(self.ref_haps)[0]
        e_mat = -np.ones((n_ref + 1, n_loci, 2))

        # Calculate Hardy-Weinberg Emissions
        e_mat[0, :, 1] = self.p  # Calculate the Emission Matrix

        # Calculate Emissions from Copying
        e_mat[1:, :, 1] = (self.ref_haps == 1)  # Copying without mistake

        # The emission of the derived Allele
        e_mat[:, :, 0] = 1 - e_mat[:, :, 1]

        # Sanity Check: Only positive E. Probabilities
        assert(np.min(e_mat) >= 0)

        if remember == True:
            self.e_mat = e_mat  # Remember it.

        return e_mat

    def give_emission_state(self, ob_stat, e_mat):
        """Gives the emission matrix of observed states
        Return emission matrix [k,l]"""
        ob_stat = ob_stat[0, :]  # Do ONLY first observed Haplotype
        assert(len(ob_stat) == np.shape(e_mat)[1])  # Sanity Check

        e_mat = e_mat[:, range(len(ob_stat)), ob_stat]
        e_mat[e_mat == 0] = self.e_rate  # Error probabilities
        e_mat[e_mat == 1] = 1 - self.e_rate
        return e_mat

    def give_emission_log(self, ob_stat):
        """Return the full emission Probability directly in Log Space.
        ob_stat: Observed Genotypes [2,l] (only use 1st row)"""
        ob_stat = ob_stat[0, :]  # Do ONLY first observed Haplotype
        assert(len(ob_stat) == np.shape(self.ref_haps)[1])  # Sanity Check

        n_loci = np.shape(self.ref_haps)[1]
        n_ref = np.shape(self.ref_haps)[0]
        e_mat0 = np.zeros((n_ref + 1, n_loci), dtype=np.float)

        # Do the HW
        e_mat0[0, :] = np.log(self.p * (1 - self.e_rate) + (1 - self.p) * self.e_rate) * (
            ob_stat == 1) + np.log((1 - self.p) * (1 - self.e_rate) + self.p * self.e_rate) * (ob_stat == 0)
        e_mat0[1:, :] = np.log(self.e_rate) * (self.ref_haps != ob_stat[None, :]) + \
            np.log(1 - self.e_rate) * (self.ref_haps == ob_stat[None, :])
        assert(np.max(e_mat0) < 0)  # Sanity Check (In Log Space Pr. <0)

        return e_mat0

###############################
###############################


class RC_Model_Emissions(Model_Emissions):
    """Implements the Read Count model Emission probabilities.
    Inherits from Model_Emission, in particular the constructor (Calculation
    of Mean Allele Frequency from the Reference and ref_haps)"""
    p = []  # Vector of mean alle frequencies in Reference [l]
    ref_haps = []  # Array of Haplotypes in Reference [n_ref, l]
    e_mat = []  # # Full Emission Matrix [n_ref+1, n_loci, 2]

    e_rate = 1e-2   # The error rate per read
    e_rate_ref = 1e-3  # The error rate for the reference genome states
    # (to not run into trouble for high coverage SNPs)

    def give_emission_matrix(self, remember=False):
        """Return Emission Matrix, which describes
        probabilities in Genotypes [n_ref+1, n_loci, 3]"""
        p = self.p
        ref_haps = self.ref_haps
        e_rate_ref = self.e_rate_ref

        n_loci = np.shape(ref_haps)[1]
        n_ref = np.shape(ref_haps)[0]
        p_hgeno = -np.ones((n_ref + 1, n_loci, 3))

        # Do the HW State 0
        p_hgeno[0, :, 0] = (1 - p) ** 2
        p_hgeno[0, :, 1] = 2 * p * (1 - p)
        p_hgeno[0, :, 2] = p ** 2

        # Do the copying states (add some error)
        p_hgeno[1:, :, 1] = e_rate_ref / 2
        p_hgeno[1:, :, 0] = (ref_haps == 0) * (1 - e_rate_ref) + \
            (ref_haps == 1) * e_rate_ref / 2
        p_hgeno[1:, :, 2] = (ref_haps == 1) * (1 - e_rate_ref) + \
            (ref_haps == 0) * e_rate_ref / 2

        # Sanity Check if genotype probabilities sum up to (approx.) 1
        assert(np.all(np.isclose(np.sum(p_hgeno, axis=2), 1)))
        assert((np.min(p_hgeno) >= 0) & (
            np.max(p_hgeno) <= 1))   # Sanity Check

        if remember == True:
            self.e_mat = p_hgeno
        return p_hgeno

    def give_emission_state(self, ob_stat, e_mat):
        """Gives the emission matrix of observed states
        Return emission matrix [n_ref+1, n_loci] of each
        ob_stat: [2, n_loci] Matrix with Nr Ref/Alt Reads in Row0/Row1 (!)
        e_mat: Probabilities of genotypes [n_ref+1, n_loci, 3]"""

        e_rate = self.e_rate  # Load the error rate per read

        # What's the probability of observing a dervided read given hidden genotypes 00 01 11
        p_read = np.array([e_rate, 0.5, 1 - e_rate])

        # Calculate the Binomial Likelihoods of RC Data
        rc_tot = np.sum(ob_stat, axis=0)
        rc_der = ob_stat[1, :]

        prob_binom = binom.pmf(
            rc_der[:, None], rc_tot[:, None], p_read[None, :])

        # Sum over each of the 3 possible genotypes
        p_full = np.sum(e_mat * prob_binom[None, :, :], axis=2)
        return p_full

    def give_emission_log(self, ob_stat):
        """Return the full emission Probability directly in Log Space.
        ob_stat: Observed Readcounts [2,l] array of 0/1 """
        e_mat = self.give_emission_matrix()
        e_mat = np.log(self.give_emission_state(ob_stat=ob_stat, e_mat=e_mat))
        assert(np.max(e_mat) < 0)  # In LOG Space (Assume Error Model)
        return e_mat

###############################
###############################


class Diploid_GT_Emissions(RC_Model_Emissions):
    """Implements the Emission probabilities for Diploid Genotype calls.
    Inherits top level Model_Emission, in particular the constructor (Calculation
    of Mean Allele Frequency from the Reference and ref_haps) and from
    RC_Model_Emission, in particular the calculation of probabilities of the
    [n_ref+1, n_loci, 3] genotype probability matrix give_emission_matrix"""
    p = []  # Vector of mean alle frequencies in Reference [l]
    ref_haps = []  # Array of Haplotypes in Reference [n_ref, l]
    e_mat = []  # # Full Emission Matrix [n_ref+1, n_loci, 2]

    e_rate = 1e-3
    e_rate_ref = 0.0  # The error rate for the reference genome states

    def give_emission_state(self, ob_stat, e_mat):
        """Gives the emission matrix of observed states
        Return emission matrix [n_ref+1, n_loci] of each
        ob_stat: [2, n_loci] Matrix with 0 or 1 for Ref or Alt
        e_mat: Probabilities of genotypes [n_ref+1, n_loci, 3]"""

        e_rate = self.e_rate  # Load the error rate per read
        nr_loci = np.shape(e_mat)[1]

        derived_gt = np.sum(ob_stat, axis=0)  # The Nr of of derived GT (0/1/2)
        assert(len(derived_gt) == nr_loci)

        # Calculate probability of ob. Genotypes Given 3 lat. genotypes [nr_loci, 3]
        # Initialise with "background error"
        prob_gt = np.ones((nr_loci, 3)) * e_rate / 2
        prob_gt[range(nr_loci), derived_gt] = 1 - e_rate

        # Sum over each of the 3 possible genotypes
        p_full = np.sum(e_mat * prob_gt[None, :, :], axis=2)
        return p_full

    def give_emission_log(self, ob_stat, dtype=np.float):
        """Return the full emission Probability directly in Log Space.
        ob_stat: Observed Readcounts [2,l] array of 0/1 """
        ref_haps = self.ref_haps
        derived_gt = np.sum(ob_stat, axis=0)  # The Nr of of derived GT (0/1/2)
        e_rate = self.e_rate  # Error Rate (Here FLIP error)

        n_loci = np.shape(self.ref_haps)[1]
        n_ref = np.shape(self.ref_haps)[0]
        assert(len(derived_gt) == n_loci)  # Sanity Check

        e_mat0 = np.zeros((n_ref + 1, n_loci), dtype=dtype)

        ### HW State
        e_mat0[0, :] = (1 - self.p) * (1 - self.p) * (derived_gt == 0) + \
            (1 - self.p) * self.p * 2 * (derived_gt == 1) + \
            self.p * self.p * (derived_gt == 2)
        ### Copying States
        e_mat0[1:, :] = (derived_gt == 0) * (ref_haps == 0) + \
                       (derived_gt == 2) * (ref_haps == 1)

        ### Do the bleeding from other states:
        e_mat0 = e_mat0 * (1 - e_rate) + (1 - e_mat0) * e_rate / 2.0
        e_mat0 = np.log(e_mat0)  # Go to Log Space

        return e_mat0

###############################
###############################
# Factory method


def load_emission_model(ref_states, e_model="haploid"):
    """Load the Emission Model"""
    if e_model == "haploid":
        e_obj = Model_Emissions(ref_states)
    elif e_model == "readcount":
        e_obj = RC_Model_Emissions(ref_states)
    elif e_model == "diploid_gt":
        e_obj = Diploid_GT_Emissions(ref_states)
    else:
        raise NotImplementedError("Emission Model not found!")
    return e_obj


#################################
# Do some testing with explicit values

if __name__ == "__main__":
    # ob_stat = np.array([[1, 5], [3, 3], [0, 2], [1, 0], [1, 1]]).T  ### RC Data
    ob_stat = np.array([[1, 1], [0, 0], [1, 0], [1, 0], [0, 1]]).T  # Diploid
    print(ob_stat)

    ref_haps = np.array([[1, 1, 1, 1, 1], [0, 0, 0, 0, 0], [1, 1, 1, 0, 0]])
    print(ref_haps)

    e_obj = load_emission_model(ref_haps, e_model="diploid_gt")  # readcount
    e_mat = e_obj.give_emission_matrix()
    e_prob = e_obj.give_emission_state(ob_stat=ob_stat, e_mat=e_mat)
    print(e_prob)
    e_prob0 = e_obj.give_emission_log(ob_stat=ob_stat)
    #e_prob1 = e_obj.give_emission_log_temp(ob_stat=ob_stat)

    print("Comparison:")
    #print(np.exp(e_prob1) - e_prob)
    #print(np.exp(e_prob1) - np.exp(e_prob0))

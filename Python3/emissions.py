"""
Class for calculating Emission Probabilities.
Contains Sub-Classes, as well as factory Method.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import numpy as np

###############################
###############################

class Emissions(object):
    """Class for emission probabilities
    Has methods to return emission probabilities"""

    def give_emission_matrix(self):
        """Return Emission Matrix"""
        raise NotImplementedError("Implement This in specific subclass.")

    def give_emission_state(self, ob_stat):
        """Gives the emission matrix of path of states"""
        raise NotImplementedError("Implement This in specific subclass.")

class Model_Emissions(Emissions):
    """Implements the haploid model Emission probabilities"""
    p = []  # Vector of alle frequencies
    ref_haps = []
    e_mat = []  # Full Emission Matrix

    e_rate = 0.01  # The Probability of an error

    def __init__(self, ref_haps=[]):
        """Initialize Class"""
        if len(ref_haps) > 0:
            self.ref_haps = ref_haps

        # Calculate the allele frequencies
        self.p = np.mean(self.ref_haps, axis=0)

    def give_emission_matrix(self, remember=True):
        """Return full Emission Matrix"""
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
        assert(len(ob_stat) == np.shape(e_mat)[1])  # Sanity Check

        e_prob = e_mat[:, range(len(ob_stat)), ob_stat]
        e_prob[e_prob == 0] = self.e_rate  # Tiny probability of emission
        e_prob[e_prob == 1] = 1 - self.e_rate
        return e_prob


##################################
#### Factory method

def load_emission_model(ref_states, e_model="haploid"):
    """Load the Emission Model"""
    if e_model == "haploid":
        e_obj = Model_Emissions(ref_states)
    elif e_model == "diploid":
        e_obj = 0  # Implement this.
    else:
        raise NotImplementedError("Emission Model not found!")

    return e_obj

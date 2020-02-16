"""
Class for calculating Transition Probabilities, i.e. 
infitesimal transition Matrices.
Contains Sub-Classes, as well as factory Method.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""
import numpy as np


class Transitions(object):
    """Class for transition probabilities.
    Has methods to return them"""
    n_ref = 0         # The Nr of reference Samples
    trans_mat = []    # The full transition Matrix
    r_map = []        # Legacy Variable

    def __init__(self, n_ref=20, r_map=[]):
        """Initialize Class"""
        self.n_ref = n_ref  # Set the Number of Reference Samples
        self.r_map = r_map  # Set the Recombination Map
        self.calc_transitions()

    def give_transitions(self):
        """Return Transition Matrix"""
        raise NotImplementedError("Implement This!")

    def calc_transitions(self, n=0):
        """Return Transition Matrix"""
        raise NotImplementedError("Implement This!")

    def set_params(self, **kwargs):
        """Set the Values."""
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.calc_transitions()  # Calculate the new transition Matrix


class Model_Transitions(Transitions):
    """Implements the Model Transitions"""
    roh_in = 0.0005     # The rate of jumping to another Haplotype
    roh_out = 0.001     # The rate of jumping out
    roh_jump = 0.02     # The rate of jumping within ROH

    def calc_transitions(self, n=0, rate=True, submat33=True):
        """Return Transition Matrix to exponate.
        n: Nr of Reference Haplotypes
        submat33: Whether to only fill in """
        if n == 0:     # Default to Nr of References as set in Class
            n = self.n_ref

        if submat33 == True:
            # Initialize Transition Matrix Only do 3 States (bc symmetry)
            t_mat = -np.ones((3, 3))
        else:
            t_mat = -np.ones((n + 1, n + 1))  # Initialize Transition Matrix

        t_mat[1:, 0] = self.roh_out  # The rate of jumping out roh
        t_mat[0, 1:] = self.roh_in / n  # Jumping into any ROH State
        t_mat[1:, 1:] = self.roh_jump / n  # Jumping between ROH State

        # Do the Diagonal (do the usual model - for inf. substract 1)
        di = np.diag_indices(np.shape(t_mat)[0])
        d = 1 - rate  # Whether to use 1 one the diagonal
        t_mat[di] = d - self.roh_out - self.roh_jump + \
            self.roh_jump / (n)  # Don't forget the self jump
        t_mat[0, 0] = d - self.roh_in   # The rate of staying in diploid

        # Sanity Check if everything was filled correctly
        if (submat33 == False) and (rate == True):
            assert(np.all(np.sum(t_mat, axis=1) > -0.0001))
            assert(np.all(np.sum(t_mat, axis=1) < 0.0001))
            #print(np.sum(t_mat, axis=1))  # Output

        self.trans_mat = t_mat
        return t_mat

    def give_transitions(self):
        """Give the transition_matrix"""
        return self.trans_mat


##################################
# Factory method

def load_transition_model(t_model="model", n_ref=20):
    """Load the Transition Model"""
    if t_model == "model":
        t_obj = Model_Transitions(n_ref=n_ref)
        # self.t_obj.set_params(n_ref=self.n_ref)

    else:
        raise NotImplementedError("Transition Model not found!")

    return t_obj

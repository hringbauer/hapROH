import numpy as np
import matplotlib.pyplot as plt
import os  # For Saving to Folder


#################################
#################################

class HMM_Analyze(object):
    """Analyze Class for HMMs.
    This is the main Class, all specific Inference schemes inherit from it
    and overwrite functions.
    Contains Parameters"""
    folder = "" # The working folder
    l = 1000  # Nr of the Observations
    n_ref = 20  # The Size of the Reference Panel
    ref_states = []  # Ref. Array of k Reference States to Copy from. [kxl]
    ob_stat = []  # The observed State [l]

    v_path = []  # The inferred Viterbi Path

    t_obj, e_obj = 0, 0  # Objects for the transitions and Emission probabilities

    def __init__(self, folder="./Simulated/Example0/",
                 t_model="model", e_model="haploid"):
        """Initialize Class"""
        self.load_data(folder=folder)
        self.load_transition_model(t_model=t_model)
        self.load_emission_model(e_model=e_model)
        self.folder = folder # Save the working folder

    def load_data(self, folder="./Simulated/Example0/"):
        """Loads the Data"""
        self.ref_states = np.loadtxt(
            folder + "refs.csv", dtype="int", delimiter=",")
        self.ob_stat = np.loadtxt(
            folder + "hap.csv", dtype="int", delimiter=",")

        # Do some Post-Processing of the Loading
        self.n_ref = np.shape(self.ref_states)[0]

        print(f"Successfully loaded Data from: {folder}")

    def load_transition_model(self, t_model="model"):
        """Load the Transition Model"""
        if t_model == "model":
            self.t_obj = Model_Transitions()

        print(f"Loaded Transition Model: {t_model}")

    def load_emission_model(self, e_model="haploid"):
        """Load the Emission Model"""
        if e_model == "haploid":
            self.e_obj = Model_Emissions(ref_haps=self.ref_states)

        if e_model == "diploid":
            self.e_obj = 0

        print(f"Loaded Emission Model: {e_model}")

    ###############################
    # The actual Inference Part

    def calc_viterbi_path(self, save=True):
        """Calculate the Viterbi Path.
        save: Boolean whether to save the Viterbi Path"""

        # 1) Get the emission and transition probabilities.
        e_mat = self.e_obj.give_emission_matrix()
        t_mat = self.t_obj.give_transitions()
        ob_stat = self.ob_stat


        e_mat0 = np.log(e_mat)  # Do the Transformation to Log Space
        t_mat0 = np.log(t_mat)

        print("Loaded Transition and Emission Matrix:")
        print(np.shape(t_mat))
        print(np.shape(e_mat))
        print("Loaded Observations:")
        print(np.shape(ob_stat))

        e_prob0 = e_mat0[:,range(len(ob_stat)), ob_stat] # The observing probabilities of the States

        n_states = np.shape(t_mat)[0]
        n_loci = np.shape(t_mat)[1]

        # Do the actual optimization (with back-tracking)
        # Initialize
        mp = np.zeros((n_states, n_loci), dtype="float")
        pt = np.zeros((n_states, n_loci), dtype="int")  # Previous State Pinter

        mp[:, 0] = np.log(0.01)  # The initial states
        mp[0, 0] = np.log(1)    # Start with no-ROH state

        for i in range(1, n_loci):  # Do the Viterbi-Iteration
            for j in range(n_states):
                new_p = mp[:, i - 1] + t_mat0[:, j] + e_prob0[j,i]
                m = np.argmax(new_p)  # Find the Maximum Probability
                mp[j, i] = new_p[m]
                pt[j, i] = m  # Set the pointer to previous path

        # Do the trace back
        path[-1] = np.argmax(mp[:, -1])  # The highest probability

        path = np.zeros(n_loci, dtype="int")
        for i in range(n_loci - 1, 0):
            path[i - 1] = pt[path[i]]  # Always th pointer to the previous path

        self.vpath = path

        print("Finished Calculation Viterbi Path:")
        print(path)
        print(f"Total log Probability: {mp[path[-1]]}")

        ## Save the Path
        if save==True:
            np.savetxt(folder + "viterbi_path.csv", path,
                       delimiter=",",  fmt='%i')  # Save the Viterbi Path

        ###############################


class Transitions(object):
    """Class for transition probabilities.
    Has methods to return them"""

    def __init__(self):
        """Initialize Class"""
        self.calc_transitions()

    def give_transitions(self):
        """Return Transition Matrix"""
        raise NotImplementedError("Implement This!")

    def calc_transitions(self, n=0):
        """Return Transition Matrix"""
        raise NotImplementedError("Implement This!")


class Model_Transitions(Transitions):
    """Implements the Model Transitions"""
    roh_in = 0.002      # The rate of jumping to another Haplotype
    roh_out = 0.01     # The rate of jumping out
    roh_jump = 0.05    # The rate of jumping within ROH
    n_ref = 20         # The Nr of reference Samples

    trans_mat = []    # The basic transition Matrix

    def calc_transitions(self, n=0):
        """Return Transition Matrix to exponate.
        n: Nr of Reference Haplotypes"""
        if n == 0:
            n = self.n_ref

        t_mat = -np.ones((n + 1, n + 1))  # Initialize Transition Matrix

        t_mat[1:, 0] = self.roh_out  # The rate of jumping out roh
        t_mat[0, 1:] = self.roh_in / n  # Jumping into any ROH State
        t_mat[1:, 1:] = self.roh_jump / n  # Jumping between ROH State

        di = np.diag_indices(n + 1)
        t_mat[di] = 1 - self.roh_out - self.roh_jump + \
            self.roh_jump / (n)  # Don't forget the self jump
        t_mat[0, 0] = 1 - self.roh_in   # The rate of staying in diploid

        # Sanity Check if everything was filled
        assert np.min(t_mat) >= 0
        # Sanity check if jumping out rates sum to 1
        assert(np.all(np.sum(t_mat, axis=1) > 0.9999))
        assert(np.all(np.sum(t_mat, axis=1) < 1.0001))
        self.trans_mat = t_mat    # Set the Transition Matrix

    def give_transitions(self):
        """Give the transition_matrix"""
        return self.trans_mat

###############################


class Emissions(object):
    """Class for emission probabilities
    Has methods to return emission probabilities"""

    def give_emission_matrix(self):
        """Return Emission Matrix"""
        raise NotImplementedError("Implement This!")


class Model_Emissions(Emissions):
    """Implements the haploid model Emission probabilities"""
    p = []  # Vector of alle frequencies
    ref_haps = []
    e_mat = []  # Full Emission Matrix

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

        assert(np.min(e_mat) >= 0)  # Assert all
        if remember == True:
            self.e_mat = e_mat  # Remember it.

        return e_mat

    def give_emission(self, pos=0, state=0):
        """Give Emission Probability"""
        raise NotImplementedError("Implement This!")


###############################
###############################
# Do some Testing
hmm = HMM_Analyze()
# hmm.load_data(folder="./Simulated/Example0/")
print(np.shape(hmm.e_obj.ref_haps))
# print(hmm.e_obj.p)
# hmm.e_obj.give_emission_matrix()
# print(np.shape(hmm.e_obj.e_mat))
# print(hmm.e_obj.e_mat[0,:,1])

hmm.calc_viterbi_path()

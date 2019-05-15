"""
Main Inference Class for HMM. Wrapper for Inerence of Posterior.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""


import numpy as np
import matplotlib.pyplot as plt
import os                     # For Saving to Folder
import cProfile               # For Profiling
from scipy.special import logsumexp
# from func import fwd_bkwd    # Import the Python Function
from cfunc import fwd_bkwd, viterbi_path, fwd_bkwd_fast  # The Cython Functions
from func import fwd_bkwd_p, viterbi_path_p   # The Python Functions
from emissions import load_emission_model     # Factory Method
from transitions import load_transition_model
from postprocessing import give_Postprocessing
from preprocessing import load_preprocessing

#################################
#################################


class HMM_Analyze(object):
    """Analyze Class for HMMs.
    This is the main Class, all specific Inference schemes inherit from it
    and overwrite functions.
    Contains Parameters"""
    folder = ""  # The working folder
    output = True
    save = True  # Whether to save output data to disk
    save_fp = True  # Whether to save the full Posterior
    n_ref = 20  # The Size of the Reference Panel [k-1]
    ref_states = []  # Ref. Array of k Reference States to Copy from. [kxl]
    ob_stat = []     # The observed State [l]

    r_map = []     # The Map position of every marker [l]

    v_path = []  # The inferred Viterbi Path [l]
    posterior = []  # The inferred Posterior Matrix [kxl] Log Space

    t_obj, e_obj = 0, 0  # Objects for Transition & Emission probabilities
    p_obj = 0          # Object that does preprocessing

    iid = ""  # Remember the Individual
    ch = 0    # Which Chromosome

    fwd_bkwd, viterbi_path = 0, 0  # Function for the fwd-bkwd Algorithm

    e_model = "haploid"
    t_model = "model"
    p_model = "SardHDF5"

    def __init__(self, folder="./Simulated/Example0/",
                 t_model="model", e_model="haploid", p_model="SardHDF5",
                 output=True, save=True, cython=True, manual_load=False,
                 save_fp=True):
        """Initialize Class. output: Boolean whether to print
        Cython: Whether to use Cython.
        Manual_Load: Whether to skip automatic loading of data"""
        self.t_model = t_model
        self.e_model = e_model
        self.p_model = p_model
        self.folder = folder  # Save the working folder
        self.output = output
        self.save = save
        self.save_fp = save_fp

        if manual_load == False:
            self.load_objects()

        # Load the Heavy Lifting Functions
        if cython == True:
            self.fwd_bkwd = fwd_bkwd
            self.viterbi_path = viterbi_path

        elif cython == 2:   # To test the "fast algorithm". Remove later
            if self.output == True:
                print("Using Linear-State Speed-Up")

            self.fwd_bkwd = fwd_bkwd_fast
            self.viterbi_path = viterbi_path

        else:
            self.fwd_bkwd = fwd_bkwd_p
            self.viterbi_path = viterbi_path_p

    def load_objects(self, iid="", ch=0, n_ref=503):
        """Load all the required Objects in right order"""
        self.load_preprocessing_model()
        self.load_data(iid, ch, n_ref)
        self.load_emission_model()
        self.load_transition_model()

    def load_data(self, iid="", ch=0, n_ref=503):
        """Load the External Data"""
        gts_ind, gts, r_map, out_folder = self.p_obj.load_data(
            iid=iid, ch=ch, n_ref=n_ref, folder=self.folder)

        self.ch = ch
        self.iid = iid
        self.ref_states = gts
        self.r_map = r_map
        self.ob_stat = gts_ind
        self.folder = out_folder

        # Do some Post-Processing for summary Parameters
        self.n_ref = np.shape(self.ref_states)[0]

        # Sanity Checks for loading
        assert(len(self.r_map) == np.shape(self.ob_stat)[1])  # Sanity Check
        assert(len(self.r_map) == np.shape(self.ref_states)[1])

        if self.output:
            print(f"Successfully loaded Data from: {self.folder}")

    def load_emission_model(self):
        """Method to load an Emission Model"""
        self.e_obj = load_emission_model(self.ref_states, e_model=self.e_model)

        if self.output:
            print(f"Loaded Emission Model: {self.e_model}")

    def load_transition_model(self):
        """Load the Transition Model"""
        self.t_obj = load_transition_model(self.t_model, n_ref=self.n_ref)

        if self.output:
            print(f"Loaded Transition Model: {self.t_model}")

    def load_preprocessing_model(self):
        self.p_obj = load_preprocessing(
            p_model=self.p_model, save=self.save, output=self.output)

        if self.output:
            print(f"Loaded Pre Processing Model: {self.p_model}")

    def set_diploid_observations(self):
        """Simulate random Diploid Observations"""
        obs = self.ob_stat
        n_loci = np.shape(obs)[1]
        phases = np.random.randint(2, size=n_loci)
        obs[0, :] = obs[phases, np.arange(n_loci)]
        self.ob_stat = obs

    def prepare_rmap(self, cm=False, min_gap=1e-10):
        """Return the recombination map
        Input: Map Positions [l]
        Return: Rec. Distance Array [l]
        cm: Whether input is in centimorgan or morgan
        min_cap: Minimum Gap between Loci"""
        r_map = np.zeros(len(self.r_map))  # Make new array
        r_map[1:] = self.r_map[1:] - self.r_map[:-1]  # Calculate Differences
        assert(np.min(r_map) >= 0)
        if cm == True:
            r_map = r_map / 100     # Normalize to CentiMorgan.
        # Extend the minimum gap where needed
        r_map = np.maximum(r_map, min_gap)

        if self.output == True:
            print(f"Minimum Genetic Map: {np.min(self.r_map):.4f}")
            print(f"Maximum Genetic Map: {np.max(self.r_map):.4f}")
            print(f"Gaps bigger than 0.1 cM: {np.sum(r_map > 0.001)}")
            print(f"Maximum Gap: {np.max(r_map) * 100:.4f} cM")

        return r_map

    def pre_compute_transition_matrix(self, t, r_vec, n_ref):
        """Precompute and return the full transition Matrix
        t full Transition Matrix [k,k]. NO LOG STATE
        r_vec Map Length of Jumps [l] """
        # n = np.shape(t)[0] - 1  # Nr of Reference States
        t_simple = prep_3x3matrix(t, n_ref)
        t_mat = exponentiate_r(t_simple, r_vec)

        # Normalize to transition rate for non-collapsed state
        # print(np.sum(t_mat, axis=2))   # Should be one: Sanity Check!

        t_mat[:, :2, 2] = t_mat[:, :2, 2] / (n_ref - 1)
        # t_mat[:, 2, 2] = t_mat[:, 1, 1]  # By Symmetr, not needed

        return t_mat

    ###############################
    ###############################
    # The actual Inference Part

    def calc_viterbi_path(self, save=True):
        """Calculate the Viterbi Path.
        save: Boolean whether to save the Viterbi Path"""

        # 1) Get the emission and transition probabilities.
        e_mat = self.e_obj.give_emission_matrix()
        t_mat = self.t_obj.give_transitions()
        # t_mat = np.eye(len(t_mat)) + t_mat LEGACY
        # Precompute the 3x3 Transition Matrix
        r_map = self.prepare_rmap()  # Get the Recombination Map
        t_mat_full = self.pre_compute_transition_matrix(
            t_mat, r_map, self.n_ref)  # [l, k, k]

        ob_stat = self.ob_stat[0, :]  # Do the first observed Haplotype

        e_prob = self.e_obj.give_emission_state(ob_stat=ob_stat, e_mat=e_mat)
        e_prob[e_prob == 0] = 1e-20  # Tiny probability of emission

        e_prob0 = np.log(e_prob)
        t_mat_full0 = np.log(t_mat_full)

        end_p = np.empty(np.shape(e_prob)[0], dtype=np.float)
        end_p[1:] = 0.001       # Low Probability
        end_p[0] = 1 - np.sum(end_p[1:])
        end_p0 = np.log(end_p)  # Go to Log Space

        if self.output == True:
            print("Loaded Transition and Emission Matrix:")
            print(np.shape(t_mat))
            print(np.shape(e_mat))
            print("Loaded Observations:")
            print(np.shape(ob_stat))

        path = self.viterbi_path(e_prob0, t_mat_full0, end_p0)

        # Save the Path
        self.vpath = path
        if save == True:
            save_path = self.folder + "viterbi_path.csv"
            np.savetxt(save_path, path,
                       delimiter=",",  fmt='%i')  # Save the Viterbi Path
            if self.output:
                print(f"Saved to: {save_path}")

        if self.output:
            print(f"Finished Calculation Viterbi Path: {path}")

    def calc_posterior(self, save=True, full=False):
        """Calculate the poserior for each path
        FULL: Wether to return fwd, bwd as well as tot_ll (Mode for postprocessing)"""
        e_mat = self.e_obj.give_emission_matrix()
        t_mat = self.t_obj.give_transitions()
        ob_stat = self.ob_stat[0, :]  # Do the first observed Haplotype
        # e_mat[e_mat == 0] = 1e-20  # Tiny probability of emission. LEGACY
        r_map = self.prepare_rmap()  # Get the Recombination Map

        n_states = np.shape(e_mat)[0]
        n_loci = np.shape(e_mat)[1]

        if self.output:
            print("Loaded Transition and Emission Matrix:")
            print(np.shape(t_mat))
            print(np.shape(e_mat))
            print("Loaded Observations:")
            print(np.shape(ob_stat))

        # The observing probabilities of the States [k,l]
        e_prob = self.e_obj.give_emission_state(ob_stat=ob_stat, e_mat=e_mat)
        e_prob0 = np.log(e_prob)

        # Initialize  the fwd and bwd probabilities
        fwd = np.ones((n_states, n_loci), dtype="float")
        fwd[:, 0] = 1e-4  # Initial Probabilities not in HW
        fwd[0, 0] = 1 - np.sum(fwd[1:, 0])  # The initial HW prob.
        fwd0 = np.log(fwd)  # Change to log space

        bwd = np.ones((n_states, n_loci), dtype="float")
        bwd[:, -1] = 1e-4  # Initial Probabilities
        bwd[0, -1] = 1 - np.sum(bwd[1:, -1])
        bwd0 = np.log(bwd)  # Change to log space

        # Precompute the 3x3 Transition Matrix
        t_mat_full = self.pre_compute_transition_matrix(
            t_mat, r_map, self.n_ref)

        # Do the forward-backward Algorithm:
        if full == False:
            post = self.fwd_bkwd(e_prob0, t_mat, fwd0, bwd0, t_mat_full)

        elif full == True:  # If FULL Mode: Return results prematurely
            post, fwd, bwd, tot_ll = self.fwd_bkwd(
                e_prob0, t_mat, fwd0, bwd0, t_mat_full, full=True)
            return post, fwd, bwd, tot_ll

        if self.output:
            print("Finished Calculation State Posteriors")

        # Save the Data
        self.posterior = post  # Remember the Posterior

        if self.save_fp == True:
            path = self.folder + "posterior.csv"
            np.savetxt(path, post,
                       delimiter=",",  fmt='%f')
            print(f"Saved Full Posterior to to {self.folder}")

        if save == True:
            path = self.folder + "posterior0.csv"
            np.savetxt(path, post[0, :],
                       delimiter=",",  fmt='%f')
            print(f"Saved Zero State Posterior to {self.folder}.")

    def optimze_ll_transition_param(self, roh_trans_params):
        """Calculate and return the log likelihoods for Transitions Parameters
        roh_trans_params [m]"""
        m = len(roh_trans_params)

        ll_hoods = []  # The Vector for the log likelihoods

        for p in roh_trans_params:
            # Set the transition Parameters
            self.t_obj.set_params(roh_jump=p)
            _, _, _, tot_ll = self.calc_posterior(save=False, full=True)
            ll_hoods.append(tot_ll)

        return np.array(ll_hoods)

    def post_processing(self, method="standard", save=True):
        """Do the Postprocessing of ROH Blocks
        Parameters: See in Postprocessing"""
        pp = give_Postprocessing(folder=self.folder, method=method,
                                 output=self.output, save=save)
        pp.call_roh(ch=self.ch, iid=self.iid)
        return


###############################
# Some Helper functions

def prep_3x3matrix(t, n_ref):
    """Prepares the grouped 3x3 Matrix (3rd State: Everything in OTHER ROH State)"""
    # n = np.shape(t)[0] - 1  # Infer the Number of reference States
    n = n_ref
    print(f"Reference Number: {n}")
    # Initiate to -1 (for later Sanity Check if everything is filled)
    t_simple = -np.ones((3, 3))
    t_simple[:2, :2] = t[:2, :2]
    # The probability of staying when in diff. ROH State:
    # t_simple[2, 2] = np.sum(t[2, 2:])   # Legacy
    t_simple[2, 2] = -np.sum(t[2, :2])  # Minus the rates of Jumping
    # print("Target Value:")
    # print(t_simple[2, 2])
    # Jumping into 3rd state: Sum over all reference states
    t_simple[:2, 2] = t[:2, 2] * (n - 1)
    t_simple[2, :2] = t[2, :2]  # The jumping out probability is the same
    # print("Total Jumping Probabilities:")
    # print(np.sum(t_simple, axis=1))
    return t_simple


def exponentiate_r(rates, rec_v):
    """Calculates exponentiation of the rates matrix with rec_v
    rates: 2D Matrix.
    rec_v: Array of length l"""
    eva, evec = np.linalg.eig(rates)  # Do the Eigenvalue Decomposition
    # Sanity Check whether Matrix is valid rate Matrix
    assert(np.max(eva) <= 1)
    evec_r = np.linalg.inv(evec)    # Do the Inversion
    # Create vector of the exponentiated diagonals
    d = np.exp(rec_v[:, None] * eva)
    # Use some Einstein Sum Convention Fun (C Speed):
    res = np.einsum('...ik, ...k, ...kj ->...ij', evec, d, evec_r)
    # Make sure that all transition rates are valuable
    assert(0 <= np.min(res))
    return res

####################################
####################################

# Joint Run for a Sardinian Sample


def analyze_individual(iid, ch, n_ref=503, save=True, save_fp=False):
    """Run the analysis for one individual and chromosome"""
    hmm = HMM_Analyze(cython=2, p_model="SardHDF5",
                      manual_load=True, save=save, save_fp=save_fp)

    hmm.load_objects(iid=iid, ch=ch, n_ref=n_ref)
    hmm.set_diploid_observations()
    hmm.t_obj.set_params(roh_in=1, roh_out=10, roh_jump=100)
    hmm.calc_viterbi_path(save=save)           # Calculate the Viterbi Path.
    hmm.calc_posterior(save=save)              # Calculate the Posterior.
    hmm.post_processing(save=save)             # Do the Post-Processing.
    print(f"Analysis of {iid} and Chr. {ch} successfully concluded!")


def profiling_run():
    """Short Function for profiling"""
    hmm = HMM_Analyze(folder="./Empirical/MA89_chr3_1000G_ROH/", cython=2)
    hmm.t_obj.set_params(roh_in=1, roh_out=10, roh_jump=100)
    # print(np.shape(hmm.e_obj.ref_haps))
    hmm.calc_posterior(save=False)


if __name__ == "__main__":
    # folder = "./Simulated/Test20r/"          # "./Simulated/Test20r/"
    # d05e e: Error Introduced. d05: Downsampled
    # folder = "./Empirical/Sard100_0-10kROH8/"
    # folder = "./Empirical/1kEUR_ROH/"
    # folder = "./Empirical/MA89_chr3_1000G_ROH/"
    folder = "./Simulated/Test20r/"
    p_model = "Folder"  # "SardHDF5" or "Folder"

    # hmm = HMM_Analyze(folder=folder, cython=2, p_model=p_model)
    # hmm.set_diploid_observations()       # Set single observation per locus.
    # hmm.t_obj.set_params(roh_in=1, roh_out=10, roh_jump=100)
    # hmm.calc_viterbi_path(save=True)           # Calculate the Viterbi Path.
    # hmm.calc_posterior(save=True)              # Calculate the Posterior.
    # hmm.post_processing(save=True)             # Do the Post-Processing.

    # ch_list = [3]   # range(1, 23)
    #ch_list = range(1, 23)

    # for ch in ch_list:
    #    print(f"Doing Chromosome: {ch}")
    #    analyze_individual(iid="I1563", ch=ch, save=True)

    analyze_individual(iid="MA89", ch=4, save=True, save_fp=True)
    #analyze_individual(iid="SEC006", ch=3, save=True, save_fp=True)

# cProfile.run('profiling_run()')
# -24816.477
# LL with correct Linkage Map: -24512.563

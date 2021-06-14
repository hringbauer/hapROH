"""
Main Inference Class for HMM. Wrapper for Inerence of Posterior.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import numpy as np
import matplotlib.pyplot as plt
import os                     # For Saving to Folder
#import psutil                 # For Memory Profiling
#import cProfile               # For Profiling
# from func import fwd_bkwd    # Import the Python Function
from hapsburg.cfunc import fwd_bkwd_fast, fwd_bkwd_lowmem, fwd_bkwd_scaled, fwd_bkwd_scaled_lowmem  # Cython Functions
from hapsburg.func import fwd_bkwd_p, sloppyROH_cumsum  # Python Functions
from hapsburg.emissions import load_emission_model     # Factory Methods
from hapsburg.transitions import load_transition_model
from hapsburg.postprocessing import load_Postprocessing
from hapsburg.preprocessing import load_preprocessing

#################################
#################################


class HMM_Analyze(object):
    """Analyze Class for HMMs.
    This is the main Class, all specific Inference schemes inherit from it
    and overwrite functions.
    Contains objects as field for pre-processing, transition as well as
    emission probabilities.
    Contains most Parameters (but some of them like the output folders are decided
    by pre-processing subclass)"""
    folder = ""  # The working folder
    output = True
    save = True  # Whether to save output data to disk
    save_fp = True  # Whether to save the full Posterior
    n_ref = 20  # The Size of the Reference Panel [k-1]
    sanity_checks = True # Can turn off for better performance. Not recomm.
    ref_states = []  # Ref. Array of k Reference States to Copy from. [kxl]
    ob_stat = []     # The observed State [l]

    r_map = []      # The Map position of every marker [l]
    pos = []        # The pysical position of every marker [l]
    v_path = []     # The inferred Viterbi Path [l]
    posterior = []  # The inferred Posterior Matrix [kxl] Log Space

    t_obj, e_obj = 0, 0  # Objects for Transition & Emission probabilities
    p_obj, post_obj = 0, 0  # Objects for pre and post-processing

    iid = ""  # Remember the Individual
    ch = 0    # Which Chromosome

    fwd_bkwd = 0  # Function for the fwd-bkwd Algorithm

    e_model = "haploid"
    t_model = "model"
    p_model = "SardHDF5"
    post_model = "Standard"

    def __init__(self, folder="./Simulated/Example0/",
                 t_model="model", e_model="haploid", p_model="SardHDF5", post_model="Standard",
                 output=True, save=True, cython=True, manual_load=False,
                 save_fp=True):
        """Initialize Class. output: Boolean whether to print
        Cython: Whether to use Cython.
        Manual_Load: Whether to skip automatic loading of data"""
        self.t_model = t_model
        self.e_model = e_model
        self.p_model = p_model
        self.post_model = post_model
        self.folder = folder  # Save the working folder
        self.output = output
        self.save = save
        self.save_fp = save_fp

        if manual_load == False:
            self.load_objects()

        # Load the Heavy Lifting Functions
        if cython == True:
            self.fwd_bkwd = fwd_bkwd_fast
            if self.output == True:
                print("Using Linear-State Cython Speed-Up")

        elif cython == 2:   # To test the "fast algorithm". Remove later
            if self.output == True:
                print("Using Low-Mem Cython Linear Speed Up.")
            self.fwd_bkwd = fwd_bkwd_lowmem #fwd_bkwd_lowmem
            
        elif cython == 3:   # To test the "fast algorithm". Remove later
            if self.output == True:
                print("Using Rescaled HMM.")
            self.fwd_bkwd = fwd_bkwd_scaled_lowmem #fwd_bkwd_scaled

        else:
            self.fwd_bkwd = fwd_bkwd_p

    def load_objects(self, iid="", ch=0):
        """Load all the required Objects in right order"""
        self.load_preprocessing_model()
        self.load_data(iid, ch)
        self.load_secondary_objects()

    def load_secondary_objects(self):
        """Load all secondary objects
        (but not the pre-processing one)"""
        self.load_emission_model()
        self.load_transition_model()
        self.load_postprocessing_model()

    def load_data(self, iid="", ch=0):
        """Load the External Data"""
        gts_ind, gts, r_map, pos, out_folder = self.p_obj.load_data(
            iid=iid, ch=ch)

        self.ch = ch
        self.iid = iid
        self.ref_states = gts
        self.r_map = r_map
        self.pos = pos
        self.ob_stat = gts_ind
        self.folder = out_folder

        ### Do some Post-Processing for summary Parameters
        self.n_ref = np.shape(self.ref_states)[0]

        ### Sanity Checks for loaded data
        if self.sanity_checks:
            assert(len(self.r_map) == np.shape(self.ob_stat)[1])  # Sanity Check
            assert(len(self.r_map) == np.shape(self.ref_states)[1])
            assert(np.min(gts_ind)>=0) # No missing genotypes
            assert(np.min(gts)>=0)

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

    def load_postprocessing_model(self):
        self.post_obj = load_Postprocessing(folder=self.folder, method=self.post_model,
                                            output=self.output, save=self.save)
        if self.output:
            print(f"Loaded Post Processing Model: {self.post_model}")

    def prepare_rmap(self, cm=False, min_gap=1e-10, max_gap=0.05):
        """Return the recombination map [in Morgan]
        Input: Map Positions [l]
        Return: Rec. Distance Array [l]
        cm: Whether input is in centimorgan or morgan
        min_cap: Minimum Gap between Loci"""
        r_map = np.zeros(len(self.r_map))  # Make new array
        r_map[1:] = self.r_map[1:] - self.r_map[:-1]  # Calculate Differences
        assert(np.min(r_map) >= 0)
        if cm == True:
            r_map = r_map / 100     # Normalize to Morgan if map in cM
        
        max_g = np.max(r_map) 
        # Extend the minimum gap where needed
        r_map = np.maximum(r_map, min_gap)
        r_map = np.minimum(r_map, max_gap)

        if self.output == True:
            print(f"Minimum Genetic Map: {np.min(self.r_map):.4f} Morgan")
            print(f"Maximum Genetic Map: {np.max(self.r_map):.4f} Morgan")
            print(f"Gaps bigger than 0.1 cM: {np.sum(r_map > 0.001)}")
            print(f"Maximum Gap: {max_g * 100:.4f} cM")
            print(f"Upper Gap Cutoff: {max_gap * 100:.4f} cM")

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

    def calc_posterior(self, save=True, full=False, in_val=1e-4):
        """Calculate the poserior for each path
        FULL: Wether to return fwd, bwd as well as tot_ll (Mode for postprocessing)
        in_val: The Initial Probability to copy from one Ind."""
        t_mat = self.t_obj.give_transitions()
        ob_stat = self.ob_stat
        r_map = self.prepare_rmap()  # Get the Recombination Map

        e_mat = self.e_obj.give_emission(ob_stat=ob_stat)
        n_states = np.shape(e_mat)[0]
        n_loci = np.shape(e_mat)[1]

        if self.output:
            print("Loaded Transition and Emission Matrix:")
            print(np.shape(t_mat))
            print(np.shape(e_mat))
            print("Loaded Observations:")
            print(np.shape(ob_stat))

        # Precompute the 3x3 Transition Matrix
        t_mat_full = self.pre_compute_transition_matrix(
            t_mat, r_map, self.n_ref)

        # Do the forward-backward Algorithm:
        if full == False:
            post = self.fwd_bkwd(e_mat, t_mat_full, in_val, full=False, output=self.output)

        elif full:  # If FULL Mode: Return results prematurely
            post, fwd, bwd, tot_ll = self.fwd_bkwd(
                e_mat, t_mat_full, in_val, full=True, output=self.output)
            return post, fwd, bwd, tot_ll

        if self.output:
            print("Finished Calculation State Posteriors")

        # Save the Data
        # self.posterior = post  # Remember the Posterior

        if self.save_fp:
            path = self.folder + "posterior.csv"
            np.savetxt(path, post,
                       delimiter=",",  fmt='%f')
            #path = self.folder + "refs.csv"
            #np.savetxt(path, self.ref_states,
            #           delimiter=",",  fmt='%i')
            print(f"Saved Full Posterior and Ref GTS to folder {self.folder}")

        if save:
            path = self.folder + "posterior0.csv"
            np.savetxt(path, post[0, :],
                       delimiter=",",  fmt='%f')
            print(f"Saved Zero State Posterior to folder {self.folder}.")

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

    def post_processing(self, save=True):
        """Do the Postprocessing of ROH Blocks
        Parameters: See in Postprocessing"""
        self.post_obj.call_roh(ch=self.ch, iid=self.iid)
        return

    def mmr_call(self, windowSize=0.001, save=True):
        """Calculate Maximal Match Rate"""
        ob_stat = self.ob_stat
        r_map = self.r_map
        ref_states = self.ref_states  # Get Reference

        max_agree_rate_window = sloppyROH_cumsum(
            r_map, ob_stat, ref_states, windowSize=windowSize)

        if save == True:
            path = self.folder + "posterior0.csv"
            np.savetxt(path, max_agree_rate_window,
                       delimiter=",",  fmt='%f')
            print(f"Saved Zero State Posterior to {path}.")

###############################
###############################
# Some Helper functions


def prep_3x3matrix(t, n_ref):
    """Prepares the grouped 3x3 Matrix (3rd State: Everything in OTHER ROH State)"""
    n = n_ref
    print(f"Reference Number: {n}")
    # Initiate to -1 (for later Sanity Check if everything is filled)
    t_simple = -np.ones((3, 3))
    t_simple[:2, :2] = t[:2, :2]
    # The probability of staying when in diff. ROH State:
    t_simple[2, 2] = -np.sum(t[2, :2])  # Minus the rates of Jumping
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


def print_memory_usage():
    """Print the current Memory Usage in mB"""
    process = psutil.Process(os.getpid())
    mb_usage = process.memory_info().rss / 1e6
    print(f"Memory Usage: {mb_usage} mB")
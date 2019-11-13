import numpy

# python sucks
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "Python3"))
# this is getting really bad
# how is python good again?
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../PackagesSupport/loadEigenstrat/"))

from hmm_inference import HMM_Analyze
from msTransitions import rate_matrix_oneSL, exponentiate_rate_matrices, allele_frequencies, one_step_transitions_reference, non_roh_state_map, non_roh_transition_matrices
# from msHMM import fancy_initial_distribution, stationary_inital_distribution, get_ROH_posterior
from msEmissions import extended_genotype_emissions, extended_binom_emission
from msFun import extended_fwd_bkwd_fast


def fake_transition_matrices (rates, rec_v):
    transition_matrices = numpy.zeros((len(rec_v),6,6))
    for i in range(len(rec_v)):
        transition_matrices[i] = numpy.eye(6)
#         transition_matrices[i] = numpy.zeros((6,6))
#         transition_matrices[i,:2,:2] = 0.5
#         transition_matrices[i,2:4,2:4] = 0.5
#         transition_matrices[i,2:4,2:4] = numpy.ones((2,2)) - numpy.eye(2)
#         transition_matrices[i,4:6,4:6] = 0.5
    return transition_matrices


def fancy_initial_distribution (n_ref, f_marg = [0.5, 0.5], pi_s = 1e-4, pi_l = 1e-4):
	"""
	Return a vector with a suitable initial distribution
	n_ref: # of reference haplotpyes
	pi_s: total weight of short ROH
	pi_l: total weight of long ROH
	"""

	init_d = numpy.zeros (4 + 2*n_ref)

	# first the non-ROH states
	pi_n = 1 - pi_l - pi_s
	(idx_to_geno, geno_to_idx) = non_roh_state_map ()
	for x in range(len(idx_to_geno)):
		init_d[x] = pi_n * f_marg[idx_to_geno[x][0]] * f_marg[idx_to_geno[x][1]]

	# then the short ROH states
	for x in range(4, 4 + n_ref):
		init_d[x] = pi_s/n_ref

	# and the long ROH states
	for x in range(4 + n_ref, 4 + 2*n_ref):
		init_d[x] = pi_l/n_ref

	return init_d


def stationary_inital_distribution (n_ref, f_marg, Q):
	"""
	Returns the initial disitrubion constructed from the stationary distribution of Q
	"""

	assert (Q.shape == (6,6))

	# first get the eigendecomposition of Q^T
	eva, evec = numpy.linalg.eig(numpy.transpose(Q))

	sortedIdx = numpy.argsort(eva)
	firstIdx = sortedIdx[-1]
	assert (numpy.isclose (eva[firstIdx], 0, rtol=1e-12, atol=1e-12))
	# and also, that there is no other 0 eigenvalue
	assert (numpy.max(eva[sortedIdx[:-1]]) < -1e-14)

	# now the first eigentvector gives the stationaray distribution
	stat = evec[:,firstIdx]/numpy.sum(evec[:,firstIdx])
	assert (numpy.min(stat) > -1e-14)

	# extract stationary mass for short ROH and long ROH
	pi_s = stat[2] + stat[3]
	pi_l = stat[4] + stat[5]

	# now just use the regular one
	return fancy_initial_distribution (n_ref, f_marg, pi_s=pi_s, pi_l=pi_l)


def get_ROH_posterior (log_full_posterior):
    """
    Takes a full posterior and sums over the ROH states to just get the ROH posterior
    """

    # exponentiate
    full_posterior = numpy.exp(log_full_posterior)

    n_loci = full_posterior.shape[0]
    n_states = full_posterior.shape[1]
    n_ref = (n_states - 4)//2

    newPost = numpy.zeros ((n_loci, 3))

    # non ROH
    newPost[:,0] = numpy.sum (full_posterior[:,0:4], axis=1)

    # short ROH
    newPost[:,1] = numpy.sum (full_posterior[:,4:(4+n_ref)], axis=1)

    # long ROH
    newPost[:,2] = numpy.sum (full_posterior[:,(4+n_ref):(4+2*n_ref)], axis=1)

    return newPost


class HapsburgFiftyThree:

    
    def __init__(self, refHaps, recoMap, in_S = 100, in_L = 1, out_S = 400, out_L = 10,
                     roh_jump = 300, e_rate_ref=1e-3, e_rate=1e-2):

        self.n_ref = refHaps.shape[0]
        self.n_loci = refHaps.shape[1]
        # need to remember this for later
        self.e_rate = e_rate

        # prepare the recombination map
        # fake class instance to use a member function of class HMM_Analyze
        fake = type('', (), {})()
        fake.r_map = recoMap
        fake.output = False
        self.r_map = HMM_Analyze.prepare_rmap (self=fake)

        # transition rate matrix Q
        # t_mat = new_calc_transitions (roh_in=100, roh_out=100, roh_jump=300)
        # self.Q = rate_matrix_oneSL (in_S = 100, in_L = 2, out_S = 400, out_L = 10, roh_jump = 300)
        self.Q = rate_matrix_oneSL (in_S = in_S, in_L = in_L, out_S = out_S, out_L = out_L, roh_jump = roh_jump)

        # marginal properbilities
        self.f_marg = allele_frequencies (refHaps)

        # initial distribution
        # self.init_d = numpy.log (fancy_initial_distribution (n_ref, f_marg[0], pi_s=1e-4, pi_l=1e-4))
        self.init_d = numpy.log (stationary_inital_distribution (self.n_ref, self.f_marg[0], self.Q))
    
        # transition
        self.transition_matrices = exponentiate_rate_matrices (rates=self.Q, rec_v=self.r_map)
        # also need the allele frequencies
        self.f_trans = one_step_transitions_reference (refHaps, self.f_marg)
        # get transition matrices between the non-ROH states
        self.non_roh_transitions = non_roh_transition_matrices (self.transition_matrices, self.f_trans, self.f_marg)

        # emission
        # prepare the emissions for all loci
        # at least as much as we can do without the actual target
        self.e_mat_geno = extended_genotype_emissions (refHaps, e_rate_ref=e_rate_ref)

        
    def compute_full_log_posterior (self, target):
        
        # now do the target specific emissions
        self.e_mat_full = numpy.log (extended_binom_emission (ob_stat=target, e_mat=self.e_mat_geno,
                                                              e_rate=self.e_rate))
        # just some checks
        assert(numpy.max(self.e_mat_full) < 0) 

        # get some posterior
        (newPost, fwd, bwd, tot_ll) = extended_fwd_bkwd_fast (init_d=self.init_d, e_prob0=self.e_mat_full,
                                          t=self.transition_matrices, f_marg=self.f_marg, f_trans=self.f_trans,
                                          non_roh_transitions=self.non_roh_transitions, full=True)

        return (newPost, fwd, bwd, tot_ll)
        
        
    def compute_reduced_posterior (self, target):
        """
        Computes the posterior for the combined states non-ROH, short-ROH, and long-ROH
        """
        
        # get posterior
        (newPost, fwd, bwd, tot_ll) = self.compute_full_log_posterior (target)
        
        # sum over the grouings of the states
        return get_ROH_posterior (newPost)

    
    def compute_long_ROH_posterior (self, target):
        """
        Computes the posterior in the long ROH state
        """
        return (self.compute_reduced_posterior (target))[:,2]
    
    
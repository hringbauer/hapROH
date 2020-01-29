import numpy as np
from scipy.special import logsumexp


def fwd_bkwd_p(e_prob0, t_mat, fwd, bwd, r_map):
    """Takes emission and transition probabilities, and calculates posteriors.
    Input: [kxl] matrices of emission, transition
    and initialized fwd and bwd probabilities. All in log Space"""
    n_states, n_loci = np.shape(e_prob0)

    # Initialize Posterior Probabilities
    post = np.zeros((n_states, n_loci), dtype="float")

    t_mat0 = np.log(np.eye(n_states) + t_mat[:,:]) # Transformation to log space

    for i in range(1, n_loci):  # Do the forward recursion
        for j in range(n_states):
            trans_ll = fwd[:, i - 1] + t_mat0[:, j]
            fwd[j, i] = e_prob0[j, i] + logsumexp(trans_ll)

    for i in range(n_loci - 1, 0, -1):  # Do the backward recursion
        for j in range(n_states):
            trans_ll = t_mat0[j, :] + e_prob0[:, i] + bwd[:, i]
            bwd[j, i - 1] = logsumexp(trans_ll)

    tot_ll = logsumexp(fwd[:, -1] + bwd[:, -1])  # Get total log likelihood

    print(f"Total Log likelihood fwd: {tot_ll: .3f}")
    # COmbine the forward and backward calculations
    post = fwd + bwd - tot_ll  # The formulat is f*b/tot_l

    return post


def viterbi_path_p(e_prob0, t_mat0, end_p0):
    """Implementation of a Viterbi Path.
    e_prob0 and t_mat0 [k,l] Matrices with Emission and Transition Probabilities.
    end_p: probability to begin/end in states [k]"""
    n_states = np.shape(e_prob0)[0]
    n_loci = np.shape(e_prob0)[1]

    # Do the actual optimization (with back-tracking)
    # Initialize
    mp = np.zeros((n_states, n_loci), dtype="float")
    pt = np.zeros((n_states, n_loci), dtype="int")  # Previous State Pinter

    mp[:,0] = end_p0

    for i in range(1, n_loci):  # Do the Viterbi-Iteration
        for j in range(n_states):
            new_p = mp[:, i - 1] + t_mat0[i, :, j] + e_prob0[j, i]
            m = np.argmax(new_p)  # Find the Maximum Probability
            mp[j, i] = new_p[m]
            pt[j, i] = m          # Set the pointer to previous path

    # Do the trace back
    path = -np.ones(n_loci, dtype="int")  # Initialize
    path[-1] = np.argmax(mp[:, -1])  # The highest probability

    for i in range(n_loci - 1, 0, -1):
        # Always th pointer to the previous path
        path[i - 1] = pt[path[i], i]

    print(f"Log likelihood Path: {mp[path[-1],-1]:.3f}")
    assert(np.min(path)>=0) #Sanity check if everything was filled up
    return path

########################################
##### The MMR Functions

def sloppyROH_cumsum(recoMap, target, refHaps, windowSize = 0.001):
    """
    Figure out per snp agreement between references and reads from target
    Return Max. fraction of reads agreeing with Ref Panel
    recoMap: Recombination Map [in Morgan] [l]
    target: Target Readcounts: [2, l]
    refHaps: Reference Haplotypes: [k, l]
    windowsize: Extension to each side [in M]
    """
    # Calculate total Nr of Mismatches for all ref
    refAgree = (1 - refHaps) * target[0, :]
    altAgree = refHaps * target[1, :]
    totalAgree = refAgree + altAgree
    
    ### what are the bounds of the window around the focal SNP
    firstNeighbor = np.searchsorted(recoMap + windowSize, recoMap)
    lastNeighbor = np.searchsorted(recoMap - windowSize, recoMap)
    lastNeighbor = np.clip(lastNeighbor, a_min=None, a_max=len(recoMap)-1) # To avoid the overshooting the last Index
    
    countsPerSnp = np.sum(target, axis=0)       # Calculate total counts per snp for normalization
    
    totalAgree_cum = np.cumsum(totalAgree, axis=1)
    pad_zeros = np.zeros((np.shape(totalAgree)[0], 1)) # To add 0s at beginning of first axis
    totalAgree_cum = np.concatenate([pad_zeros, totalAgree_cum], axis=1) 
    countsPerSnp_cum = np.r_[0, np.cumsum(countsPerSnp)] # Add 0 in beginning
    
    totalAgree_window = totalAgree_cum[:, lastNeighbor] - totalAgree_cum[:, firstNeighbor]
    countsPerSnp_window = countsPerSnp_cum[lastNeighbor] - countsPerSnp_cum[firstNeighbor]
    
    agree_rate_window = totalAgree_window / countsPerSnp_window    # Average Agree Rate per Window
    max_agree_rate_window = np.max(agree_rate_window, axis=0)   # For the "best" Ref Haplotype
    return max_agree_rate_window












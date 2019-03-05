import numpy as np
from scipy.special import logsumexp


def fwd_bkwd(e_prob0, t_mat0, fwd, bwd):
    """Takes emission and transition probabilities, and calculates posteriors.
    Input: [kxl] matrices of emission, transition
    and initialized fwd and bwd probabilities. All in log Space"""
    n_states, n_loci = np.shape(e_prob0)

    # Initialize Posterior Probabilities
    post = np.zeros((n_states, n_loci), dtype="float")

    for i in range(1, n_loci):  # Do the forward recursion
        for j in range(n_states):
            trans_ll = fwd[:, i - 1] + t_mat0[:, j]
            fwd[j, i] = e_prob0[j, i] + logsumexp(trans_ll)

    for i in range(n_loci - 1, 0, -1):  # Do the backward recursion
        for j in range(n_states):
            trans_ll = t_mat0[j, :] + e_prob0[:, i] + bwd[:, i]
            bwd[j, i - 1] = logsumexp(trans_ll)

    tot_ll = logsumexp(fwd[:, -1] + bwd[:, -1])  # Get total log likelihood
    #tot_llb = logsumexp(fwd[:, 0] + bwd[:, 0])

    print(f"Total Log likelihood fwd: {tot_ll: .3f}")
    # print(f"Total Log likelihood bwd: {tot_llb: .3f})

    # COmbine the forward and backward calculations
    post = fwd + bwd - tot_ll  # The formulat is f*b/tot_l

    return post

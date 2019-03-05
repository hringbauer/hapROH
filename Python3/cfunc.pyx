# cython: language_level=3, boundscheck=False
import numpy as np
cimport cython

from scipy.special import logsumexp

DTYPE = np.float # The float data type
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.

def fwd_bkwd(double[:, :] e_prob0, double[:, :] t_mat0,
    double[:, :] fwd, double[:, :] bwd):
    """Takes emission and transition probabilities, and calculates posteriors.
    Input: kxl matrices of emission, transition
    and initialized fwd and bwd probabilities. All in log Space"""
    cdef int n_states = e_prob0.shape[0]
    cdef int n_loci = e_prob0.shape[1]
    cdef Py_ssize_t i, j, k # The Array Indices

    # Initialize Posterior and Transition Probabilities
    post = np.empty([n_states, n_loci], dtype=DTYPE)

    trans_ll = np.empty(n_states, dtype=DTYPE)
    cdef double[:] trans_ll_view = trans_ll

    for i in range(n_loci-1):  # Do the forward recursion
        i+=1
        for j in range(n_states):
          for k in range(n_states):
            trans_ll_view[k] = fwd[k, i - 1] + t_mat0[k, j]
          fwd[j, i] = e_prob0[j, i] + logsumexp(trans_ll_view)

    for i in range(n_loci-1, 0, -1):  # Do the backward recursion
      for j in range(n_states):
        for k in range(n_states):
          trans_ll_view[k] = t_mat0[j, k] + e_prob0[k, i] + bwd[k, i]
        bwd[j, i - 1] = logsumexp(trans_ll_view)

    # Get total log likelihood
    tot_ll = logsumexp(np.asarray(fwd[:, -1]) + np.asarray(bwd[:, -1]))
    #tot_llb = logsumexp(fwd[:, 0] + bwd[:, 0])

    print(f"Total Log likelihood fwd: {tot_ll: .3f}")
    # print(f"Total Log likelihood bwd: {tot_llb: .3f})

    # Combine the forward and backward calculations
    post = np.asarray(fwd) + np.asarray(bwd)
    post = post - tot_ll
    return post

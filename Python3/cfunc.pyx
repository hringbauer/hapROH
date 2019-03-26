# cython: language_level=3, boundscheck=False, wraparound=False
import numpy as np
cimport cython

#from scipy.special import logsumexp
from libc.math cimport exp, log   # For doing Exponentials and Logs

DTYPE = np.float # The float data type

cdef inline double logsumexp(double[:] vec):
  """Do the Log of the Sum of Exponentials."""
  cdef Py_ssize_t i  # The iterator Variable
  cdef double result = 0.0
  cdef double largest = vec[0]

  for i in range(1, vec.shape[0]):   # Find Maximum in vec
      if (vec[i] > largest):
          largest = vec[i]
  for i in range(vec.shape[0]):
      result += exp(vec[i] - largest)
  return largest + log(result)


def fwd_bkwd(double[:, :] e_prob0, double[:, :] t_mat0,
    double[:, :] fwd, double[:, :] bwd, full=False):
    """Takes emission and transition probabilities, and calculates posteriors.
    Input: kxl matrices of emission, transition
    and initialized fwd and bwd probabilities. Given in log Space
    full: Boolean whether to return everything"""
    cdef int n_states = e_prob0.shape[0]
    cdef int n_loci = e_prob0.shape[1]
    cdef Py_ssize_t i, j, k # The Array Indices


    # Initialize Posterior and Transition Probabilities
    post = np.empty([n_states, n_loci], dtype=DTYPE)

    trans_ll = np.empty(n_states, dtype=DTYPE)
    cdef double[:] trans_ll_view = trans_ll

    for i in range(1, n_loci):  # Do the forward recursion
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
    for k in range(n_states):  # Simply sum the two 1D arrays
      trans_ll_view[k] = fwd[k, n_states-1] + bwd[k, n_states-1]
    tot_ll = logsumexp(trans_ll_view)

    print(f"Total Log likelihood: {tot_ll: .3f}")

    # Combine the forward and backward calculations
    fwd1 = np.asarray(fwd)  # Transform
    bwd1 = np.asarray(bwd)
    post = fwd1 + bwd1 - tot_ll

    if full==False:
      return post

    elif full==True:   # Return everything
      return post, fwd1, bwd1, tot_ll


def fwd_bkwd_fast(double[:, :] e_prob0, double[:, :] t_mat0,
    double[:, :] fwd, double[:, :] bwd, full=False):
    """Takes emission and transition probabilities, and calculates posteriors.
    Input: kxl matrices of emission, transition
    and initialized fwd and bwd probabilities. Given in log Space
    full: Boolean whether to return everything. Uses speed-up specific for Genotype data
    (with many of the same transition rates)"""
    cdef int n_states = e_prob0.shape[0]
    cdef int n_loci = e_prob0.shape[1]
    cdef Py_ssize_t i, j, k    # The Array Indices
    cdef double stay           # The Probablility of Staying

    # Initialize Posterior and Transition Probabilities
    post = np.empty([n_states, n_loci], dtype=DTYPE)

    trans_ll = np.empty(n_states-1, dtype=DTYPE) # Array for pre-calculations
    cdef double[:] trans_ll_view = trans_ll

    three_v = np.empty(3, dtype=DTYPE)     # Array of size three
    cdef double[:] three_v_view = three_v

    two_v = np.empty(2, dtype=DTYPE)       # Array of size two
    cdef double[:] two_v_view = two_v

    # Precalculate "stay" term
    stay = np.log(np.exp(t_mat0[1, 1]) - np.exp(t_mat0[1, 2]))
    for i in range(1, n_loci):  # Do the forward recursion
        for k in range(1, n_states): # Calculate logsum of ROH states:
            trans_ll_view[k-1] = fwd[k, i - 1]
        f_l = logsumexp(trans_ll_view) # Logsum of ROH States

        # Do the 0 State:
        two_v_view[0] = fwd[0, i - 1] + t_mat0[0, 0]   # Staying in 0 State
        two_v_view[1] = f_l + t_mat0[1, 0]             # Going into 0 State
        fwd[0, i] = e_prob0[0, i] + logsumexp(two_v_view)

        ### Do the other states
        # Preprocessing:
        three_v_view[0] = fwd[0, i - 1] + t_mat0[0, 1]   # Coming from 0 State
        three_v_view[1] = f_l + t_mat0[1,2]             # Coming from other ROH State

        for j in range(1, n_states):  # Do the final run over all states
          three_v_view[2] = fwd[j, i-1] +  stay
          fwd[j, i] = e_prob0[j, i] + logsumexp(three_v_view)

    #############################
    ### Do the Backward Algorithm
    for i in range(n_loci-1, 0, -1):  # Do the backward recursion
      for k in range(1, n_states): # Calculate logsum of ROH states:
          trans_ll_view[k-1] = bwd[k, i] + e_prob0[k, i]
      f_l = logsumexp(trans_ll_view) # Logsum of ROH States

      # Do the 0 State:
      two_v_view[0] = bwd[0, i] + t_mat0[0, 0] + e_prob0[0, i]   # Staying in 0 State
      two_v_view[1] = f_l + t_mat0[0, 1]             # Going into 0 State
      bwd[0, i - 1] = logsumexp(two_v_view)

      ### Do the other states
      # Preprocessing:
      three_v_view[0] = e_prob0[0, i] + bwd[0, i] + t_mat0[1, 0]  # Coming from 0 State
      three_v_view[1] = f_l + t_mat0[1,2]             # Coming from other ROH State

      for j in range(1, n_states):  # Do the final run over all states
        three_v_view[2] = e_prob0[j, i] + bwd[j, i] +  stay
        bwd[j, i - 1] = logsumexp(three_v_view)  # Fill in the backward Probability

    # Get total log likelihood
    for k in range(n_states):  # Simply sum the two 1D arrays
      trans_ll_view[k] = fwd[k, n_states-1] + bwd[k, n_states-1]
    tot_ll = logsumexp(trans_ll_view)

    print(f"Total Log likelihood: {tot_ll: .3f}")

    # Combine the forward and backward calculations
    fwd1 = np.asarray(fwd)  # Transform
    bwd1 = np.asarray(bwd)
    post = fwd1 + bwd1 - tot_ll

    if full==False:
      return post

    elif full==True:   # Return everything
      return post, fwd1, bwd1, tot_ll


def viterbi_path(double[:, :] e_prob0, double[:, :] t_mat0, double[:] end_p0):
    """Implementation of a Viterbi Path.
    e_prob0 and t_mat0 [k,l] Matrices with Emission and Transition Probabilities.
    end_p: probability to begin/end in states [k]"""
    cdef int n_states = e_prob0.shape[0]
    cdef int n_loci = e_prob0.shape[1]
    cdef Py_ssize_t i, j, k # The Array Indices
    cdef int m  # Placeholder for Maximum
    cdef double v # Value to set

    # Do the actual optimization (with back-tracking)
    # Initialize Views:
    cdef double[:, :] mp = np.empty((n_states, n_loci), dtype=np.float)
    cdef double[:] new_p = np.empty(n_states, dtype = np.float) # Temporary Array
    cdef long[:,:] pt = np.empty((n_states, n_loci), dtype = np.int)  # Previous State Pointer

    for k in range(n_states):
      mp[k, 0] = end_p0[k]

    for i in range(1, n_loci):  # Do the Viterbi-Iteration
        for j in range(n_states):
          for k in range(n_states):
              new_p[k] = mp[k, i - 1] + t_mat0[k, j] + e_prob0[j, i]

          m, v = 0, new_p[0]    # Own Implementation of Argmax
          for k in range(1,n_states):   # Find Maximum
            if (new_p[k] > v):
              m, v = k, new_p[k]

          mp[j, i] = new_p[m]
          pt[j, i] = m          # Set the pointer to previous path

    # Do the trace back
    cdef long[:] path = -np.ones(n_loci, dtype=np.int)  # Initialize

    x = np.argmax(mp[:, n_loci-1])  # The highest probability
    path[n_loci-1] = x

    for i in range(n_loci - 1, 0, -1):
        # Always th pointer to the previous path
        path[i - 1] = pt[path[i], i]

    m = path[n_loci-1]
    print(f"Log likelihood Path: {mp[m,n_loci-1]:.3f}")
    assert(np.min(path)>=0) #Sanity check if everything was filled up
    return np.asarray(path)

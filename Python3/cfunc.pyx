# cython: language_level=3, boundscheck=True, wraparound=False
import numpy as np
import psutil      # For Memory Profiling
import os
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

cdef inline long argmax(double[:] vec):
  """Return Max and ArgMax"""
  cdef Py_ssize_t i  # The iterator Variable.
  cdef int m = 0     # Position of the Maximum.
  cdef double v = vec[0]

  for k in range(1, vec.shape[0]):   # Find Maximum
    if (vec[k] > v):
      m, v = k, vec[k]
  return m  # Return Argmax


def print_memory_usage():
    """Print the current Memory Usage in mB"""
    process = psutil.Process(os.getpid())
    mb_usage = process.memory_info().rss / 1e6
    print(f"Memory Usage: {mb_usage} mB")


def fwd_bkwd(double[:, :] e_prob0, double[:, :] t_mat,
    double[:, :] fwd, double[:, :] bwd, double[:,:,:] t, full=False):
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

    ### Transform to Log space
    cdef double[:, :] t_mat0 = np.log(np.eye(n_states) + t_mat[:,:])  # Do log of (relevant) transition Matrix

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
      trans_ll_view[k] = fwd[k, n_loci-1] + bwd[k, n_loci-1]
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


def fwd_bkwd_fast(double[:, :] e_prob0, double[:, :] t_mat,
    double[:, :, :] t, double in_val = 1e-4, full=False):
    """Takes emission and transition probabilities, and calculates posteriors.
    Uses speed-up specific for Genotype data (pooling same transition rates)
    Input:
    Emission probabilities [k x l] (log space)       (log space)
    Transition probabilities (infinitesimal) [k x k] (normal space)
    Initialized fwd and bwd probabilities [k x l]    (log space)
    t: Transition Matrix: [l x 3 x 3]                     (normal space)
    full: Boolean whether to return everything"""
    cdef int n_states = e_prob0.shape[0]
    cdef int n_loci = e_prob0.shape[1]
    cdef Py_ssize_t i, j, k    # The Array Indices
    cdef double stay           # The Probablility of Staying

    # Initialize Posterior and Transition Probabilities
    post = np.empty([n_states, n_loci], dtype=DTYPE)

    trans_ll = np.empty(n_states-1, dtype=DTYPE) # Array for pre-calculations
    cdef double[:] trans_ll_view = trans_ll

    trans_ll1 = np.empty(n_states, dtype=DTYPE) # Array for calculations
    cdef double[:] trans_ll_view1 = trans_ll1

    three_v = np.empty(3, dtype=DTYPE)     # Array of size three
    cdef double[:] three_v_view = three_v

    two_v = np.empty(2, dtype=DTYPE)       # Array of size two
    cdef double[:] two_v_view = two_v

    # Do transform to Log Space:
    cdef double[:,:,:] t0 = np.log(t)         # Do log of recombination Map

    ### Initialize FWD BWD matrices
    fwd0 = np.zeros((n_states, n_loci), dtype="float")
    fwd0[:, 0] = np.log(in_val)  # Initial Probabilities
    fwd0[0, 0] = np.log(1 - (n_states - 1) * in_val)
    cdef double[:,:] fwd = fwd0

    bwd0 = np.zeros((n_states, n_loci), dtype="float")
    bwd0[:, -1] = np.log(in_val)
    bwd0[0, -1] = np.log(1 - (n_states - 1) * in_val)
    cdef double[:,:] bwd = bwd0

    #############################
    ### Do the Forward Algorithm
    for i in range(1, n_loci):  # Run forward recursion
        stay = log(t[i, 1, 1] - t[i, 1, 2])  # Do the log of the Stay term

        for k in range(1, n_states): # Calculate logsum of ROH states:
            trans_ll_view[k-1] = fwd[k, i - 1]
        f_l = logsumexp(trans_ll_view) # Logsum of ROH States

        # Do the 0 State:
        two_v_view[0] = fwd[0, i - 1] + t0[i, 0, 0]   # Staying in 0 State
        two_v_view[1] = f_l + t0[i, 1, 0]             # Going into 0 State
        fwd[0, i] = e_prob0[0, i] + logsumexp(two_v_view)

        ### Do the other states
        # Preprocessing:
        three_v_view[0] = fwd[0, i - 1] + t0[i, 0, 1]   # Coming from 0 State
        three_v_view[1] = f_l + t0[i, 1, 2]             # Coming from other ROH State

        for j in range(1, n_states):  # Do the final run over all states
          three_v_view[2] = fwd[j, i-1] +  stay
          fwd[j, i] = e_prob0[j, i] + logsumexp(three_v_view)

    #############################
    ### Do the Backward Algorithm
    for i in range(n_loci-1, 0, -1):  # Run backward recursion
      stay = log(t[i, 1, 1] - t[i, 1, 2])

      for k in range(1, n_states): # Calculate logsum of ROH states:
          trans_ll_view[k-1] = bwd[k, i] + e_prob0[k, i]
      f_l = logsumexp(trans_ll_view) # Logsum of ROH States

      # Do the 0 State:
      two_v_view[0] = bwd[0, i] + t0[i, 0, 0] + e_prob0[0, i]   # Staying in 0 State
      two_v_view[1] = f_l + t0[i, 0, 1]                         # Going into 0 State
      bwd[0, i - 1] = logsumexp(two_v_view)

      ### Do the other states
      # Preprocessing:
      three_v_view[0] = e_prob0[0, i] + bwd[0, i] + t0[i, 1, 0]
      three_v_view[1] = f_l + t0[i, 1, 2]    # Coming from other ROH State

      for j in range(1, n_states):  # Do the final run over all states
        three_v_view[2] = e_prob0[j, i] + bwd[j, i] +  stay
        bwd[j, i - 1] = logsumexp(three_v_view)  # Fill in the backward Probability

    # Get total log likelihood
    for k in range(n_states):  # Simply sum the two 1D arrays
      trans_ll_view1[k] = fwd[k, n_loci - 1] + bwd[k, n_loci - 1]
    tot_ll = logsumexp(trans_ll_view1)

    # Combine the forward and backward calculations
    fwd1 = np.asarray(fwd, dtype=np.float)  # Transform
    bwd1 = np.asarray(bwd, dtype=np.float)
    post = fwd1 + bwd1 - np.float(tot_ll)

    print("Memory Usage Full:")
    print_memory_usage()   ## For MEMORY_BENCH
    print(f"Total Log likelihood: {tot_ll: .3f}")

    if full==False:
      return post

    elif full==True:   # Return everything
      return post, fwd1, bwd1, tot_ll


def fwd_bkwd_lowmem(double[:, :] e_prob0, double[:, :] t_mat,
    double[:, :, :] t, double in_val = 1e-4, full=False):
    """Takes emission and transition probabilities, and calculates posteriors.
    Uses speed-up specific for Genotype data (pooling same transition rates)
    Low-Mem: Do no save the full FWD BWD and Posterior. Use temporary
    Arrays for saving.
    Input:
    Emission probabilities [k x l] (log space)       (log space)
    Transition probabilities (infinitesimal) [k x k] (normal space)
    Initialized fwd and bwd probabilities [k x l]    (log space)
    t: Transition Matrix: [l x 3 x 3]                     (normal space)
    full: Boolean whether to return everything"""
    cdef int n_states = e_prob0.shape[0]
    cdef int n_loci = e_prob0.shape[1]
    cdef Py_ssize_t i, j, k    # The Array Indices
    cdef double stay           # The Probablility of Staying
    cdef double tot_ll  # The total Likelihood (need for Posterior)

    # Initialize Posterior and Transition Probabilities
    post = np.empty(n_loci, dtype=DTYPE) # Array of 0 State Posterior
    cdef double[:] post_view = post

    trans_ll = np.empty(n_states-1, dtype=DTYPE) # Array for pre-calculations
    cdef double[:] trans_ll_view = trans_ll

    trans_ll1 = np.empty(n_states, dtype=DTYPE) # Array for calculations
    cdef double[:] trans_ll_view1 = trans_ll1

    three_v = np.empty(3, dtype=DTYPE)     # Array of size three
    cdef double[:] three_v_view = three_v

    two_v = np.empty(2, dtype=DTYPE)       # Array of size two
    cdef double[:] two_v_view = two_v

    # Do transform to Log Space:
    cdef double[:,:,:] t0 = np.log(t)         # Do log of recombination Map

    ### Initialize FWD BWD Arrays
    fwd0 = np.zeros(n_states, dtype=DTYPE)
    fwd0[:] = np.log(in_val)  # Initial Probabilities
    fwd0[0] = np.log(1 - (n_states - 1) * in_val)
    cdef double[:] fwd = fwd0

    bwd0 = np.zeros(n_states, dtype=DTYPE)
    bwd0[:] = np.log(in_val)
    bwd0[0] = np.log(1 - (n_states - 1) * in_val)
    cdef double[:] bwd = bwd0

    tmp0 = np.zeros(n_states, dtype=DTYPE)
    cdef double[:] tmp = tmp0

    #############################
    ### Do the Forward Algorithm

    post_view[0] = fwd[0] # Add to first locus 0 Posterior
    for i in range(1, n_loci):  # Run forward recursion
        stay = log(t[i, 1, 1] - t[i, 1, 2])  # Do the log of the Stay term

        for k in range(1, n_states): # Calculate logsum of ROH states:
            trans_ll_view[k-1] = fwd[k]
        f_l = logsumexp(trans_ll_view) # Logsum of ROH States

        # Do the 0 State:
        two_v_view[0] = fwd[0] + t0[i, 0, 0]   # Staying in 0 State
        two_v_view[1] = f_l + t0[i, 1, 0]             # Going into 0 State
        tmp[0] = e_prob0[0, i] + logsumexp(two_v_view)

        ### Do the other states
        # Preprocessing:
        three_v_view[0] = fwd[0] + t0[i, 0, 1]   # Coming from 0 State
        three_v_view[1] = f_l + t0[i, 1, 2]             # Coming from other ROH State

        for j in range(1, n_states):  # Do the final run over all states
          three_v_view[2] = fwd[j] +  stay
          tmp[j] = e_prob0[j, i] + logsumexp(three_v_view)

        ### Make tmp new fwd vec:
        for j in range(0, n_states):
          fwd[j] = tmp[j]
        post_view[i] = fwd[0]  # Add to 0-State Posterior

    ### Get total log likelihood
    for k in range(n_states):  # Simply sum the two 1D arrays
      trans_ll_view1[k] = fwd[k] + bwd[k]
    tot_ll = logsumexp(trans_ll_view1)

    #############################
    ### Do the Backward Algorithm
    ## last0-State Posterior
    post_view[n_loci-1] = post_view[n_loci-1] + bwd[0] - tot_ll

    for i in range(n_loci-1, 0, -1):  # Run backward recursion
      stay = log(t[i, 1, 1] - t[i, 1, 2])

      for k in range(1, n_states): # Calculate logsum of ROH states:
          trans_ll_view[k-1] = bwd[k] + e_prob0[k, i]
      f_l = logsumexp(trans_ll_view) # Logsum of ROH States

      # Do the 0 State:
      two_v_view[0] = bwd[0] + t0[i, 0, 0] + e_prob0[0, i]   # Staying in 0 State
      two_v_view[1] = f_l + t0[i, 0, 1]                         # Going into 0 State
      tmp[0] = logsumexp(two_v_view)

      ### Do the other states
      # Preprocessing:
      three_v_view[0] = e_prob0[0, i] + bwd[0] + t0[i, 1, 0]
      three_v_view[1] = f_l + t0[i, 1, 2]    # Coming from other ROH State

      for j in range(1, n_states):  # Do the final run over all states
        three_v_view[2] = e_prob0[j, i] + bwd[j] +  stay
        tmp[j] = logsumexp(three_v_view)  # Fill in the backward Probability

      ### Make tmp new bwd vec:
      for j in range(0, n_states):
        bwd[j] = tmp[j]

      ### Finalize the 0 Posterior
      post_view[i-1] = post_view[i-1] + bwd[0] - tot_ll

    print(f"Total Log likelihood: {tot_ll: .3f}")
    print_memory_usage()   ## For MEMORY_BENCH

    if full==False:
      return post[None,:]  # For "fake" axis

    elif full==True:   # Return everything
      return post[None,:], fwd0, bwd0, tot_ll


def viterbi_path(double[:, :] e_prob0, double[:, :, :] t_mat0, double[:] end_p0):
    """Implementation of a Viterbi Path.
    e_prob0  Matrices with Emission Probabilities, [k,l] (log space)
    t_mat: Transition Matrix: [l x 3 x 3]  (normal space)
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

    trans_ll = np.empty(n_states-1, dtype=DTYPE) # Array for pre-calculations
    cdef double[:] trans_ll_view = trans_ll

    three_v = np.empty(3, dtype=DTYPE)     # Array of size three
    cdef double[:] three_v_view = three_v

    three_vi = np.empty(3, dtype=int)       # Int Array of size three
    cdef long[:] three_vi_view = three_vi

    two_v = np.empty(2, dtype=DTYPE)       # Array of size two
    cdef double[:] two_v_view = two_v

    two_vi = np.empty(2, dtype=int)       # Int Array of size two
    cdef long[:] two_vi_view = two_vi

    for k in range(n_states):
      mp[k, 0] = end_p0[k]  # Initialize with Ending Probabilities

      two_vi_view[0] = 0
      three_vi_view[0] = 0

    for i in range(1, n_loci):  # Do the Viterbi-Iteration
        ### Precomputation:
        # Do the maximal log probability of 1, ...k State:
        m = argmax(mp[1:, i - 1])
        v = mp[m+1, i - 1]

        # Do the States from collapsed states
        two_vi_view[1] = m+1   # Set the Pointers
        three_vi_view[1] = m+1

        two_v_view[1] = v + t_mat0[i, 1, 0]
        three_v_view[1] = v + t_mat0[i, 1, 2] # Move in from other ROH

        ### Do the zero State
        two_v_view[0] = mp[0, i - 1] + t_mat0[i, 0, 0]

        m = argmax(two_v_view)      ### Do a Maximum
        v = two_v_view[m]
        mp[0, i] = v + e_prob0[0, i]   ### Set Max. Probability
        pt[0, i] = two_vi_view[m]      ### Set Pointer for Backtrace

        ### Do the other States
        three_v_view[0] = mp[0, i - 1] + t_mat0[i, 0, 1] # Move from 0 State

        for k in range(1, n_states):   # Find Maximum
          three_v_view[2] = mp[k, i - 1] + t_mat0[i, 1, 1] # The Stay State
          three_vi_view[2] = k

          m = argmax(three_v_view)      ### Do a Maximum
          v = three_v_view[m]

          mp[k, i] = v + e_prob0[k, i]   ### Set Max. Probability
          pt[k, i] = three_vi_view[m]      ### Set Pointer for Backtrace

    ### Do the trace back
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

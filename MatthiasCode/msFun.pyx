# cython: language_level=3, boundscheck=True, wraparound=False
import numpy
import psutil      # For Memory Profiling
import os
cimport cython

#from scipy.special import logsumexp
from libc.math cimport exp, log   # For doing Exponentials and Logs

DTYPE = numpy.float # The float data type

cdef inline double logsumexp(double[:] vec):
  """Do the Log of the Sum of Exponentials."""
  cdef Py_ssize_t i  # The iterator Variable
  cdef double result = 0.0
  cdef double largest = vec[0]

  for i in range(1, vec.shape[0]):   # Find Maximum in vec
      if (vec[i] > largest):
          largest = vec[i]
  if (largest == float("-inf")):
      return largest
  for i in range(vec.shape[0]):
      result += exp(vec[i] - largest)
  return largest + log(result)

# cdef inline double logdiff (double a, double b):
#   """Do the Log of the Sum of Exponentials."""
#   # cdef Py_ssize_t i  # The iterator Variable
#   cdef double result = 0.0
#   # cdef double largest = vec[0]

#   # catch zeros
#   if (a == float("-inf")):
#     assert (b == float("-inf"))
#     return float("-inf")
#   # just make sure its not too bad
#   assert (a > b - 1e-12)
#   if (a <= b):
#     # we get zero
#     return float("-inf")
#   # else do stuff
#   # normalize by a
#   return a + log(1 - exp(b - a))  

cdef inline long argmax(double[:] vec):
  """Return Max and ArgMax"""
  cdef Py_ssize_t i  # The iterator Variable.
  cdef int m = 0     # Position of the Maximum.
  cdef double v = vec[0]

  for k in range(1, vec.shape[0]):   # Find Maximum
    if (vec[k] > v):
      m, v = k, vec[k]
  return m  # Return Argmax

cdef inline int getFirstAllele (int gdx):
  """Return first allele of genetype with index gdx"""
  if ((gdx == 0) or (gdx == 1)):
    return 0
  else :
    return 1

cdef inline int getSecondAllele (int gdx):
  """Return second allele of genetype with index gdx"""
  if ((gdx == 0) or (gdx == 2)):
    return 0
  else :
    return 1


def print_memory_usage():
    """Print the current Memory Usage in mB"""
    process = psutil.Process(os.getpid())
    mb_usage = process.memory_info().rss / 1e6
    print(f"Memory Usage: {mb_usage} mB")


def old_fwd_bkwd_fast(double[:, :] e_prob0, double[:, :, :] t, double in_val = 1e-4, full=False):
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
    post = numpy.empty([n_states, n_loci], dtype=DTYPE)

    trans_ll = numpy.empty(n_states-1, dtype=DTYPE) # Array for pre-calculations
    cdef double[:] trans_ll_view = trans_ll

    trans_ll1 = numpy.empty(n_states, dtype=DTYPE) # Array for calculations
    cdef double[:] trans_ll_view1 = trans_ll1

    three_v = numpy.empty(3, dtype=DTYPE)     # Array of size three
    cdef double[:] three_v_view = three_v

    two_v = numpy.empty(2, dtype=DTYPE)       # Array of size two
    cdef double[:] two_v_view = two_v

    # Do transform to Log Space:
    cdef double[:,:,:] t0 = numpy.log(t)         # Do log of recombination Map

    ### Initialize FWD BWD matrices
    fwd0 = numpy.zeros((n_states, n_loci), dtype="float")
    fwd0[:, 0] = numpy.log(in_val)  # Initial Probabilities
    fwd0[0, 0] = numpy.log(1 - (n_states - 1) * in_val)
    cdef double[:,:] fwd = fwd0

    bwd0 = numpy.zeros((n_states, n_loci), dtype="float")
    bwd0[:, -1] = numpy.log(in_val)
    bwd0[0, -1] = numpy.log(1 - (n_states - 1) * in_val)
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
    fwd1 = numpy.asarray(fwd, dtype=numpy.float)  # Transform
    bwd1 = numpy.asarray(bwd, dtype=numpy.float)
    post = fwd1 + bwd1 - numpy.float(tot_ll)

    print("Memory Usage Full:")
    print_memory_usage()   ## For MEMORY_BENCH
    print(f"Total Log likelihood: {tot_ll: .3f}")

    if full==False:
      return post

    elif full==True:   # Return everything
      return post, fwd1, bwd1, tot_ll


def new_fwd_bkwd_fast(double[:, :] e_prob0, double[:, :, :] t, double in_val = 1e-4, full=False):
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
    post = numpy.empty([n_states, n_loci], dtype=DTYPE)

    trans_ll = numpy.empty(n_states-1, dtype=DTYPE) # Array for pre-calculations
    cdef double[:] trans_ll_view = trans_ll

    trans_ll1 = numpy.empty(n_states, dtype=DTYPE) # Array for calculations
    cdef double[:] trans_ll_view1 = trans_ll1

    three_v = numpy.empty(3, dtype=DTYPE)     # Array of size three
    cdef double[:] three_v_view = three_v

    two_v = numpy.empty(2, dtype=DTYPE)       # Array of size two
    cdef double[:] two_v_view = two_v

    # Do transform to Log Space:
    cdef double[:,:,:] t0 = numpy.log(t)         # Do log of recombination Map

    ### Initialize FWD BWD matrices
    fwd0 = numpy.zeros((n_states, n_loci), dtype="float")
    fwd0[:, 0] = numpy.log(in_val)  # Initial Probabilities
    fwd0[0, 0] = numpy.log(1 - (n_states - 1) * in_val)
    cdef double[:,:] fwd = fwd0

    bwd0 = numpy.zeros((n_states, n_loci), dtype="float")
    # bwd0[:, -1] = numpy.log(in_val)
    # bwd0[0, -1] = numpy.log(1 - (n_states - 1) * in_val)
    # init bwd with 0 = log(1)
    bwd0[:, -1] = 0
    bwd0[0, -1] = 0
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
    fwd1 = numpy.asarray(fwd, dtype=numpy.float)  # Transform
    bwd1 = numpy.asarray(bwd, dtype=numpy.float)
    post = fwd1 + bwd1 - numpy.float(tot_ll)

    print("Memory Usage Full:")
    print_memory_usage()   ## For MEMORY_BENCH
    print(f"Total Log likelihood: {tot_ll: .3f}")

    if full==False:
      return post

    elif full==True:   # Return everything
      return post, fwd1, bwd1, tot_ll


def extended_fwd_bkwd_fast (double[:] init_d, double[:, :] e_prob0, double[:, :, :] t, double[:, :] f_marg, double[:, :] f_trans,
                            double[:, :, :] non_roh_transitions, double pi_s = 1e-4, double pi_l = 1e-4, full=False):
    """
    Takes emission and transition probabilities, and calculates posteriors.
    Uses speed-up specific for Genotype data (pooling same transition rates)
    Input:
    e_prob0: Emission probabilities [k x l] (log space)       (log space)
    t: Transition Matrix: [l x 6 x 6]                     (normal space)
    f_marg: Marginal allele frequencies [l x 2]
    f_trans: One step allele transitions [l x 4]
    non_roh_transitions: Transition probabilities btween non-roh [l x 4] (normal space)
    pi_s: stationary prob in short ROH
    pi_l: stationary prob in long ROH
    full: Boolean whether to return everything
    """
    cdef int n_states = e_prob0.shape[1]
    cdef int n_loci = e_prob0.shape[0]
    cdef int n_ref = (n_states - 4)//2
    cdef int begin_non_roh = 0
    cdef int begin_short_roh = 4
    cdef int begin_long_roh = begin_short_roh + n_ref
    cdef Py_ssize_t i, j, k          # The Array Indices
    # cdef double stay_short           # The Probablility of Staying
    # cdef double stay_long            # The Probablility of Staying
    cdef double stat_geno            # stationary probability for certain genotype
    # cdef double tmp
    # cdef double swap

    # Initialize Posterior and Transition Probabilities
    post = numpy.empty([n_loci, n_states], dtype=DTYPE)

    trans_ll = numpy.empty(n_ref, dtype=DTYPE) # Array for pre-calculations
    cdef double[:] trans_ll_view = trans_ll

    # trans_ll1 = numpy.empty(n_states, dtype=DTYPE) # Array for calculations
    # cdef double[:] trans_ll_view1 = trans_ll1

    three_v = numpy.empty(3, dtype=DTYPE)     # Array of size three
    cdef double[:] three_v_view = three_v

    # two_v = numpy.empty(2, dtype=DTYPE)       # Array of size two
    # cdef double[:] two_v_view = two_v

    four_v = numpy.empty(4, dtype=DTYPE)       # Array of size six
    cdef double[:] four_v_view = four_v

    six_v = numpy.empty(6, dtype=DTYPE)       # Array of size six
    cdef double[:] six_v_view = six_v

    # make a copy of the intial distribution, cause we gonna mess with it later
    initial_dist = numpy.copy (init_d)
    cdef double[:] initial_dist_view = initial_dist
    # Transform to Log Space:
    # cdef double[:,:,:] t0 = numpy.log(t)         # Do log of recombination Map
    old_settings = numpy.seterr (divide='ignore', invalid='ignore')
    cdef double[:,:,:] log_non_roh_transitions = numpy.log (non_roh_transitions)
    # print (log_non_roh_transitions[2,3,3])
    # three_v[0] = -0.5
    # three_v[1] = -0.5
    # three_v[2] = log_non_roh_transitions[2,3,3]
    # print (logsumexp(three_v))
    numpy.seterr(**old_settings)
    # make sure it all looks good
    assert (begin_long_roh + n_ref == n_states)

    # always useful to have intial probabilities around
    # first the non-ROH
    # for y in range(0, 4):
    #   stat_geno_y = log(f_marg[0,getFirstAllele(y)]*f_marg[0,getSecondAllele(y)])
    #   initial_dist[y] =  numpy.log(1 - pi_s - pi_l) + stat_geno_y
    # # all short ROH
    # initial_dist[begin_short_roh:begin_long_roh] = numpy.log (pi_s/n_ref)
    # # all long ROH
    # initial_dist[begin_long_roh:] = numpy.log (pi_l/n_ref)
    # renormalize (cause numerics can be slightly off)
    # initial_dist[:] -= logsumexp(initial_dist)
    # print (numpy.exp(logsumexp(initial_dist)))

    ### Initialize FWD BWD matrices
    fwd0 = numpy.zeros((n_loci, n_states), dtype=DTYPE)

    # Initial Probabilities
    # stationary
    fwd0[0, :] = initial_dist
    # add (log multpliy)
    fwd0[0, :] += e_prob0[0,:]
    # nice view
    cdef double[:,:] fwd = fwd0

    bwd0 = numpy.zeros((n_loci, n_states), dtype=DTYPE)
    # all last entried zero (= log(1))
    bwd0[-1,:] = 0
    cdef double[:,:] bwd = bwd0

    #############################
    ### Do the Forward Algorithm
    for l in range(1, n_loci):  # Run forward algorithm
        # print(f"++++++++++++++++++++ Next locus {l}")
        # print(numpy.exp(logsumexp(fwd0[l-1,:])))

        # were we not in ROH?
        for k in range(0,4):
          four_v_view[k] = fwd[l - 1, k]
        non_roh_lm1 = logsumexp(four_v_view)
        # print (numpy.exp(non_roh_lm1))

        # were we in short ROH?
        for k in range(0, n_ref):
            trans_ll_view[k] = fwd[l - 1, begin_short_roh + k]
        short_lm1 = logsumexp(trans_ll_view) # Logsum of short ROH States
        # print (numpy.exp(short_lm1))

        # were we in long ROH?
        for k in range(0, n_ref):
            trans_ll_view[k] = fwd[l - 1, begin_long_roh + k]
        long_lm1 = logsumexp(trans_ll_view) # Logsum of short ROH States
        # print (numpy.exp(long_lm1))
        # print (numpy.sum(numpy.exp([non_roh_lm1, short_lm1, long_lm1])))

        # transitions for non-ROH States:
        for y in range(0, 4):
          stat_geno_y = log(f_marg[l,getFirstAllele(y)]*f_marg[l,getSecondAllele(y)])
          # print (numpy.exp(stat_geno_y))
          for x in range(0, 4):
            # coming from x
            six_v_view[x] = fwd[l-1, x] + log_non_roh_transitions[l, x, y]
            # print (six_v_view[x])
          # coming from short
          six_v_view[4] = short_lm1 + log((t[l, 2, 0] + t[l, 2, 1])) + stat_geno_y
          # print (six_v_view[4])
          # coming from long
          six_v_view[5] = long_lm1 + log((t[l, 4, 0] + t[l, 4, 1])) + stat_geno_y
          # print (six_v_view[5])
          # and put it all together
          fwd[l, y] = e_prob0[l, y] + logsumexp(six_v_view)
          # print (logsumexp(six_v_view))
          # print (fwd[l, y])

        ### short ROH states
        # Preprocessing:
        four_v_view[0] = non_roh_lm1 + log((t[l, 0, 2] + t[l, 0, 3])/n_ref)   # Coming from non-ROH
        four_v_view[1] = short_lm1 + log(t[l, 2, 3]/n_ref)                    # Coming from other short ROH State (stay later)
        four_v_view[2] = long_lm1 + log((t[l, 4, 2] + t[l, 4, 3])/n_ref)      # Coming from other long ROH State
        # print (numpy.exp(four_v_view[1]))

        # stay_short = log(t[l, 2, 2] - t[l, 2, 3]/n_ref)  # Do the log of the Stay term
        # print (numpy.exp(stay_short))
        for j in range(begin_short_roh, begin_short_roh + n_ref):  # Do the final run over all short ROH states
          # swap = four_v_view[1]
          # four_v_view[1] = logdiff (four_v_view[1], fwd[l-1, j] + log(t[l, 2, 3]/n_ref))
          # four_v_view[3] = fwd[l-1, j] +  stay_short
          four_v_view[3] = fwd[l-1, j] + log(t[l, 2, 2])
          # print (numpy.exp(fwd[l-1, j])
          # print (numpy.exp(four_v_view[1]))
          # print (numpy.exp(four_v_view[3]))
          fwd[l, j] = e_prob0[l, j] + logsumexp(four_v_view)
          # trans_ll_view[j-begin_short_roh] = fwd[l, j]
          # four_v_view[1] = swap
        # print (numpy.exp(logsumexp(trans_ll_view)))

        ### long ROH states
        # Preprocessing:
        four_v_view[0] = non_roh_lm1 + log((t[l, 0, 4] + t[l, 0, 5])/n_ref)   # Coming from non-ROH
        four_v_view[1] = short_lm1 + log((t[l, 2, 4] + t[l, 2, 5])/n_ref)     # Coming from other short ROH State
        four_v_view[2] = long_lm1 + log(t[l, 4, 5]/n_ref)                     # Coming from other long ROH State (stay later)

        # stay_long = log(t[l, 4, 4] - t[l, 4, 5]/n_ref)  # Do the log of the Stay term
        for j in range(begin_long_roh, begin_long_roh + n_ref):  # Do the final run over all long ROH states
          four_v_view[3] = fwd[l-1, j] +  log(t[l, 4, 4])
          fwd[l, j] = e_prob0[l, j] + logsumexp(four_v_view)


    #############################
    ### Do the Backward Algorithm
    for l in range(n_loci-2, -1, -1):  # Run backward recursion
      # stay = log(t[i, 1, 1] - t[i, 1, 2])
      # print(f"++++++++++++++++++++ Next locus {l}")

      # sum of non-ROH states (with emission)
      for y in range(0,4):
        stat_geno_y = log(f_marg[l+1,getFirstAllele(y)]*f_marg[l+1,getSecondAllele(y)])
        four_v_view[y] = bwd[l + 1, y] + e_prob0[l+1, y] + stat_geno_y
      non_roh_lp1 = logsumexp(four_v_view)

      # sum of short ROH states (with emission)
      for k in range(begin_short_roh, begin_short_roh +n_ref):
          trans_ll_view[k-begin_short_roh] = bwd[l + 1, k] + e_prob0[l+1, k]
      short_lp1 = logsumexp(trans_ll_view) # Logsum of short ROH States
      # print (numpy.exp(short_lm1))

      # sum of long ROH (with emission)
      for k in range(begin_long_roh, begin_long_roh + n_ref):
          trans_ll_view[k-begin_long_roh] = bwd[l + 1, k] + e_prob0[l+1, k]
      long_lp1 = logsumexp(trans_ll_view) # Logsum of short ROH States

      # transitions for non-ROH States:
      for x in range(0, 4):
        # print (numpy.exp(stat_geno_y))
        for y in range(0, 4):
          # coming from x
          six_v_view[y] = bwd[l+1, y] + log_non_roh_transitions[l+1, x, y] + e_prob0[l+1, y]
          # print (six_v_view[y])
        # coming from short
        six_v_view[4] = short_lp1 + log((t[l+1, 0, 2] + t[l+1, 0, 3])/n_ref)
        # print (six_v_view[4])
        # coming from long
        six_v_view[5] = long_lp1 + log((t[l+1, 0, 4] + t[l+1, 0, 5])/n_ref)
        # print (six_v_view[5])
        # and put it all together
        bwd[l, x] = logsumexp(six_v_view)
        # print (logsumexp(six_v_view))
        # print (fwd[l, y])

      ### short ROH states
      # Preprocessing:
      four_v_view[0] = non_roh_lp1 + log((t[l+1, 2, 0] + t[l+1, 3, 0]))       # going to non-ROH
      four_v_view[1] = short_lp1 + log(t[l+1, 2, 3]/n_ref)                    # going to other short ROH State (stay later)
      four_v_view[2] = long_lp1 + log((t[l+1, 2, 4] + t[l+1, 2, 4])/n_ref)    # going to other long ROH State
      # print (numpy.exp(four_v_view[1]))

      # stay_short = log(t[l, 2, 2] - t[l, 2, 3]/n_ref)  # Do the log of the Stay term
      # print (numpy.exp(stay_short))
      for j in range(begin_short_roh, begin_short_roh + n_ref):  # Do the final run over all short ROH states
        # swap = four_v_view[1]
        # four_v_view[1] = logdiff (four_v_view[1], fwd[l-1, j] + log(t[l, 2, 3]/n_ref))
        # four_v_view[3] = fwd[l-1, j] +  stay_short
        # stay in short ROH
        four_v_view[3] = bwd[l+1, j] + log(t[l+1, 2, 2]) + e_prob0[l+1, j]
        # print (numpy.exp(fwd[l-1, j])
        # print (numpy.exp(four_v_view[1]))
        # print (numpy.exp(four_v_view[3]))
        bwd[l, j] = logsumexp(four_v_view)
        # trans_ll_view[j-begin_short_roh] = fwd[l, j]
        # four_v_view[1] = swap
      # print (numpy.exp(logsumexp(trans_ll_view)))

      ### long ROH states
      # Preprocessing:
      four_v_view[0] = non_roh_lp1 + log((t[l+1, 4, 0] + t[l+1, 5, 0]))         # going to non-ROH
      four_v_view[1] = short_lp1 + log((t[l+1, 4, 2] + t[l+1, 4, 3])/n_ref)     # going to other short ROH State
      four_v_view[2] = long_lp1 + log(t[l+1, 4, 5]/n_ref)                       # going to other long ROH State (stay later)

      # stay_long = log(t[l, 4, 4] - t[l, 4, 5]/n_ref)  # Do the log of the Stay term
      for j in range(begin_long_roh, begin_long_roh + n_ref):  # Do the final run over all long ROH states
        # stay in long ROH
        four_v_view[3] = bwd[l+1, j] +  log(t[l+1, 4, 4]) + e_prob0[l+1, j]
        bwd[l, j] = logsumexp(four_v_view)

    # # Get total log likelihood
    # for k in range(n_states):  # Simply sum the two 1D arrays
    #   trans_ll_view1[k] = fwd[k, n_loci - 1] + bwd[k, n_loci - 1]
    # for now
    # log_likelihood should just be sum over last
    tot_ll_fwd = logsumexp(fwd0[-1,:])
    # print (numpy.exp(tot_ll_fwd))
    # also get total ll from backward
    for i in range(0, n_states):
      # ruin initial dist =)
      initial_dist_view[i] += bwd[0,i]
      initial_dist_view[i] += e_prob0[0,i]
    tot_ll_bwd = logsumexp(initial_dist_view)
    # print (numpy.exp(tot_ll_bwd))
    assert (numpy.isclose (tot_ll_fwd, tot_ll_bwd, rtol=1e-12, atol=1e-12))

    tot_ll = tot_ll_fwd
    # Combine the forward and backward calculations
    fwd1 = numpy.asarray(fwd, dtype=DTYPE)  # Transform
    bwd1 = numpy.asarray(bwd, dtype=DTYPE)
    post = fwd1 + bwd1 - numpy.float(tot_ll)

    print("Memory Usage Full:")
    print_memory_usage()   ## For MEMORY_BENCH
    # print(f"Total Log likelihood: {tot_ll: .3f}")

    if full==False:
      return post

    elif full==True:   # Return everything
      return post, fwd1, bwd1, tot_ll


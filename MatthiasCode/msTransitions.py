import numpy


def old_calc_transitions (n, roh_in=100, roh_out=100, roh_jump=300, rate=True, submat33=True):
    """Return Transition Matrix to exponate.
    n: Nr of Reference Haplotypes
    submat33: Whether to only fill in """
    # if n == 0:     # Default to Nr of References as set in Class
    #     n = self.n_ref

    if submat33 == True:
        # Initialize Transition Matrix Only do 3 States (bc symmetry)
        t_mat = -numpy.ones((3, 3))
    else:
        t_mat = -numpy.ones((n + 1, n + 1))  # Initialize Transition Matrix

    t_mat[1:, 0] = roh_out  # The rate of jumping out roh
    t_mat[0, 1:] = roh_in / n  # Jumping into any ROH State
    t_mat[1:, 1:] = roh_jump / n  # Jumping between ROH State

    # Do the Diagonal (do the usual model - for inf. substract 1)
    di = numpy.diag_indices(numpy.shape(t_mat)[0])
    d = 1 - rate  # Whether to use 1 one the diagonal
    t_mat[di] = d - roh_out - roh_jump + \
        roh_jump / (n)  # Don't forget the self jump
    t_mat[0, 0] = d - roh_in   # The rate of staying in diploid

    # Sanity Check if everything was filled correctly
    if (submat33 == False) and (rate == True):
        assert(numpy.all(numpy.sum(t_mat, axis=1) > -0.0001))
        assert(numpy.all(numpy.sum(t_mat, axis=1) < 0.0001))
        #print(numpy.sum(t_mat, axis=1))  # Output

    return t_mat

def new_calc_transitions (roh_in=100, roh_out=100, roh_jump=300):
    """Return Transition Matrix to exponate.
    n: Nr of Reference Haplotypes
    submat33: Whether to only fill in """
    # if n == 0:     # Default to Nr of References as set in Class
    #     n = self.n_ref

    t_mat = -numpy.ones((3, 3))

    t_mat[1:, 0] = roh_out  # The rate of jumping out roh
    t_mat[0, 1:] = roh_in/2  # Jumping into any ROH State
    t_mat[1:, 1:] = roh_jump  # Jumping between ROH State

    # Do the Diagonal (do the usual model - for inf. substract 1)
    di = numpy.diag_indices(numpy.shape(t_mat)[0])
    t_mat[di] = - roh_out - roh_jump
    # this is at end, cause overwritten by other stuff
    t_mat[0, 0] = - roh_in   # The rate of staying in diploid

    # print (t_mat)
    assert(numpy.all(numpy.sum(t_mat, axis=1) > -0.0001))
    assert(numpy.all(numpy.sum(t_mat, axis=1) < 0.0001))
    #print(numpy.sum(t_mat, axis=1))  # Output

    return t_mat


def rate_matrix_oneSL (in_S = 100, in_L = 2, out_S = 400, out_L = 10, roh_jump = 200):
    """
    6x6 transition matrix for model with 2 ROH states and 1-locus memory for non-ROH.
    in_S: Rate to jump into short ROH.
    in_L: Rate to jump into long ROH.
    out_S: Rate to jump out of short ROH.
    out_L: Rate to jump out of short ROH.
    roh_jump: Rate to switch copying haplotpe.
    """

    Q = -numpy.zeros ((6, 6))

    # the haplotype switching rates within blocks
    Q[0,1] = Q[1,0] = roh_jump
    Q[2,3] = Q[3,2] = roh_jump
    Q[4,5] = Q[5,4] = roh_jump

    # jumping into ROH
    # short
    Q[0,2] = Q[0,3] = Q[1,2] = Q[1,3] = 0.5 * in_S
    # long
    Q[0,4] = Q[0,5] = Q[1,4] = Q[1,5] = 0.5 * in_L

    # jumping out of ROH
    # short
    Q[2,0] = Q[2,1] = Q[3,0] = Q[3,1] = 0.5 * out_S
    # long
    Q[4,0] = Q[4,1] = Q[5,0] = Q[5,1] = 0.5 * out_L

    # and make the diagonal minus the off diagonals
    di = numpy.diag_indices (numpy.shape(Q)[0])
    Q[di] = - numpy.sum (Q, axis=1)

    return Q


def exponentiate_rate_matrices (rates, rec_v):
    """
    Calculates exponentiation of the rates matrix with rec_v
    rates: 2D Matrix.
    rec_v: Array of length l
    """
    eva, evec = numpy.linalg.eig(rates)  # Do the Eigenvalue Decomposition
    # Sanity Check whether Matrix is valid rate Matrix
    # maybe use machine revision here too
    assert(numpy.max(eva) <= 1e13)
    evec_r = numpy.linalg.inv(evec)    # Do the Inversion
    # Create vector of the exponentiated diagonals
    d = numpy.exp(rec_v[:, None] * eva)
    # Use some Einstein Sum Convention Fun (C Speed):
    res = numpy.einsum('...ik, ...k, ...kj ->...ij', evec, d, evec_r)
    # Make sure that all transition rates are valid
    # as long as nothing goes more negativ then machine precision, we are fine
    assert(-1e-15 <= numpy.min(res))
    # but we have to cut it off
    res[res <= 0] = 0
    # perhaps renormalize?
    # res[]
    return res


def allele_frequencies (refHaps):
    """
    Computes allele frequencies at each locus.
    refHaps: nxL genotype matrix for reference panel (0 = ref allele, 1 = alt allele)
    Returns: 2xL matrix with allele freq of ref [0,:] and alr [1,:]
    """

    n_loci = refHaps.shape[1]
    n_ref = refHaps.shape[0]

    f = numpy.zeros ((n_loci, 2))

    # allele freq for alternative allele
    f[:,1] = numpy.sum (refHaps, axis=0)/n_ref
    # allele freq for reference allele
    f[:,0] = 1 - f[:,1]

    return f


def non_roh_state_map ():
    idx_to_geno = [(0,0),(0,1),(1,0),(1,1)]
    geno_to_idx = {}
    for i in range(len(idx_to_geno)):
        geno_to_idx[idx_to_geno[i]] = i

    return (idx_to_geno, geno_to_idx)


def one_step_transitions_reference (refHaps, f_marg):
    """
    Computes one step transitions of alleles in the reference panel.
    """
    newShape = (refHaps.shape[0], refHaps.shape[1]-1, 2)
    twoLocusHaps = - numpy.ones (newShape)

    # copy refHaps so that we have the two locus haps in place
    twoLocusHaps[:,:,0] = refHaps[:,:-1]
    twoLocusHaps[:,:,1] = refHaps[:,1:]

    # multiply first position by two, to get the indicees right
    twoLocusHaps[:,:,0] = 2 * twoLocusHaps[:,:,0]

    # and then sum it to finally get the indicees
    twoLocusHapsIdx = numpy.sum (twoLocusHaps, axis=2)

    # now we need to know at each locus, the probability of having a certain index
    f_joint = numpy.zeros ((refHaps.shape[1], 4))
    # put in the marginal at first locus
    f_joint[0,0] = f_joint[0,2] = f_marg[0,0]
    f_joint[0,1] = f_joint[0,3] = f_marg[0,1]
    # and fill the rest
    # for loops slow =(
    for l in range(1, f_joint.shape[0]):
        (idxs, counts) = numpy.unique (twoLocusHapsIdx[:,l-1], return_counts=True)
        f_joint[l,idxs.astype(int)] = counts/numpy.sum(counts)
        
    # and now divide by marginal to get the transition probabilities
    # no warnings
    invalid_err = numpy.geterr()['invalid']
    numpy.seterr (invalid='ignore')
    f_trans = numpy.copy(f_joint)
    f_trans[1:,:2] = f_trans[1:,:2]/f_marg[:-1,0][:,None]
    f_trans[1:,2:] = f_trans[1:,2:]/f_marg[:-1,1][:,None]
    # back to warnings
    numpy.seterr (invalid=invalid_err)

    # fill in the marginal distribution where there is nan
    # first put zeros
    nan_mask = numpy.isnan(f_trans)
    f_trans[numpy.isnan(f_trans)] = 0
    # figure out new stuff
    new_stuff = numpy.zeros (f_trans.shape)
    new_stuff[:,0] =  f_marg[:,0]
    new_stuff[:,1] =  f_marg[:,1]
    new_stuff[:,2] =  f_marg[:,0]
    new_stuff[:,3] =  f_marg[:,1]
    # make zeros where we don't need it
    new_stuff[numpy.invert(nan_mask)] = 0
    # and add new stuff to f_trans
    f_trans = f_trans + new_stuff

    return f_trans


def non_roh_transition_matrices (transition_matrices, f_trans, f_marg):
    # get the map
    (idx_to_geno, geno_to_idx) = non_roh_state_map()
    assert (len(idx_to_geno) == 4)
    # here come the transition rates between the non-ROH states
    non_roh_transitions = - numpy.ones ((transition_matrices.shape[0], 4, 4))
    # first one just get's the identity
    non_roh_transitions[0,:,:] = numpy.eye(4)
    # now the rest
    # scary nested for loops
    for idx in range(4):
        # get the geno for this idx
        x = idx_to_geno[idx]
    #     print (x)
        for jdx in range(4):
            # get the geno for this jdx
            y = idx_to_geno[jdx]
    #         print(y)
            # and fill the transitions
            first_transition = geno_to_idx[(x[0],y[0])]
            second_transition = geno_to_idx[(x[1],y[1])]
    #         print ("++++", first_transition, second_transition)
            non_roh_transitions[1:,idx,jdx] = transition_matrices[1:,0,0] * f_trans[1:,first_transition] * f_trans[1:,second_transition]
            non_roh_transitions[1:,idx,jdx] += transition_matrices[1:,0,1]/2 * f_trans[1:,first_transition] * f_marg[1:,y[1]]
            non_roh_transitions[1:,idx,jdx] += transition_matrices[1:,0,1]/2 * f_marg[1:,y[0]] * f_trans[1:,second_transition]
  
    return non_roh_transitions


def HW_transition_matrices (transition_matrices, f_marg):
    # get the map
    (idx_to_geno, geno_to_idx) = non_roh_state_map()
    assert (len(idx_to_geno) == 4)
    # here come the transition rates between the non-ROH states
    non_roh_transitions = - numpy.ones ((transition_matrices.shape[0], 4, 4))
    # first one just get's the identity
    non_roh_transitions[0,:,:] = numpy.eye(4)
    # now the rest
    # scary nested for loops
    for idx in range(4):
        # get the geno for this idx
        x = idx_to_geno[idx]
    #     print (x)
        for jdx in range(4):
            # get the geno for this jdx
            y = idx_to_geno[jdx]
    #         print(y)
            # and fill the transitions
            # just HW
            non_roh_transitions[1:,idx,jdx] = (transition_matrices[1:,0,0] + transition_matrices[1:,0,1]) *f_marg[1:,y[0]] * f_marg[1:,y[1]]
  
    return non_roh_transitions


def johnsQ (fast_rate, slow_rate, switch_rate):
    Q = numpy.zeros ((6,6))
    
    # the haplotype switching rates within blocks
    Q[2,3] = Q[3,2] = fast_rate
    Q[4,5] = Q[5,4] = slow_rate

    # jumping between ROH states
    Q[2,4] = Q[2,5] = Q[3,4] = Q[3,5] = switch_rate
    Q[4,2] = Q[4,3] = Q[5,2] = Q[5,3] = switch_rate

    # and make the diagonal minus the off diagonals
    di = numpy.diag_indices (numpy.shape(Q)[0])
    Q[di] = - numpy.sum (Q, axis=1)

    return Q


    

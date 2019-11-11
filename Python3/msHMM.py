import numpy

from msTransitions import non_roh_state_map


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


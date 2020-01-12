import numpy
from scipy.stats import binom  # Binomial Likelihood
# it seems like importing from same directory works
# only need the state map
from msTransitions import non_roh_state_map

def give_emission_matrix (ref_haps, e_rate_ref):
    """Return Emission Matrix, which describes
    probabilities in Genotypes [n_ref+1, n_loci, 3]"""
    p = numpy.mean (ref_haps, axis=0)
    n_ref = numpy.shape(ref_haps)[0]
    n_loci = numpy.shape(ref_haps)[1]
    p_hgeno = -numpy.ones((n_ref + 1, n_loci, 3))

    # Do the HW State 0
    p_hgeno[0, :, 0] = (1 - p) ** 2
    p_hgeno[0, :, 1] = 2 * p * (1 - p)
    p_hgeno[0, :, 2] = p ** 2

    # Do the copying states (add some error)
    p_hgeno[1:, :, 1] = e_rate_ref / 2
    p_hgeno[1:, :, 0] = (ref_haps == 0) * (1 - e_rate_ref) + \
        (ref_haps == 1) * e_rate_ref / 2
    p_hgeno[1:, :, 2] = (ref_haps == 1) * (1 - e_rate_ref) + \
        (ref_haps == 0) * e_rate_ref / 2
        
    # Sanity Check if genotype probabilities sum up to (approx.) 1
    assert(numpy.all(numpy.isclose(numpy.sum(p_hgeno, axis=2), 1)))
    assert((numpy.min(p_hgeno) >= 0) & (
        numpy.max(p_hgeno) <= 1))   # Sanity Check

    return p_hgeno


def give_emission_state (ob_stat, e_mat, e_rate):
    """Gives the emission matrix of observed states
    Return emission matrix [n_ref+1, n_loci] of each
    ob_stat: [2, n_loci] Matrix with Nr Ref/Alt Reads in Row0/Row1 (!)
    e_mat: Probabilities of genotypes [n_ref+1, n_loci, 3]"""

    # What's the probability of observing a dervided read given hidden genotypes 00 01 11
    p_read = numpy.array([e_rate, 0.5, 1 - e_rate])

    # Calculate the Binomial Likelihoods of RC Data
    rc_tot = numpy.sum(ob_stat, axis=0)
    rc_der = ob_stat[1, :]

    prob_binom = binom.pmf(
        rc_der[:, None], rc_tot[:, None], p_read[None, :])

    # Sum over each of the 3 possible genotypes
    p_full = numpy.sum(e_mat * prob_binom[None, :, :], axis=2)
    return p_full


def old_calc_rc_emission_log (ref_haps, ob_stat, e_rate_ref=1e-3, e_rate=1e-2):
    """Return the full emission Probability directly in Log Space.
    ob_stat: Observed Readcounts [2,l] array of 0/1 """
    e_mat = give_emission_matrix (ref_haps, e_rate_ref)
    e_mat = numpy.log (give_emission_state (ob_stat=ob_stat, e_mat=e_mat, e_rate=e_rate))
    assert(numpy.max(e_mat) < 0)  # In LOG Space (Assume Error Model)
    return e_mat


def new_give_emission_matrix (ref_haps, e_rate_ref):
    """Return Emission Matrix, which describes
    probabilities in Genotypes [n_ref+1, n_loci, 3]"""
    p = numpy.mean (ref_haps, axis=0)
    n_ref = numpy.shape(ref_haps)[0]
    n_loci = numpy.shape(ref_haps)[1]
    p_hgeno = -numpy.ones((n_ref + 1, n_loci, 3))

    # Do the HW State 0
    p_hgeno[0, :, 0] = (1 - p) ** 2
    p_hgeno[0, :, 1] = 2 * p * (1 - p)
    p_hgeno[0, :, 2] = p ** 2

    # Do the copying states (add some error)
    p_hgeno[1:, :, 1] = 2*e_rate_ref*(1-e_rate_ref)
    p_hgeno[1:, :, 0] = (ref_haps == 0) * numpy.power(1 - e_rate_ref,2) + \
        (ref_haps == 1) * numpy.power(e_rate_ref,2)
    p_hgeno[1:, :, 2] = (ref_haps == 1) * numpy.power(1 - e_rate_ref,2) + \
        (ref_haps == 0) * numpy.power(e_rate_ref,2)
        
    # Sanity Check if genotype probabilities sum up to (approx.) 1
    assert(numpy.all(numpy.isclose(numpy.sum(p_hgeno, axis=2), 1)))
    assert((numpy.min(p_hgeno) >= 0) & (
        numpy.max(p_hgeno) <= 1))   # Sanity Check

    return p_hgeno


def new_calc_rc_emission_log (ref_haps, ob_stat, e_rate_ref=1e-3, e_rate=1e-2):
    """Return the full emission Probability directly in Log Space.
    ob_stat: Observed Readcounts [2,l] array of 0/1 """
    e_mat = new_give_emission_matrix (ref_haps, e_rate_ref)
    e_mat = numpy.log (give_emission_state (ob_stat=ob_stat, e_mat=e_mat, e_rate=e_rate))
    assert(numpy.max(e_mat) < 0)  # In LOG Space (Assume Error Model)
    return e_mat


def extended_rc_emission_log (ref_haps, ob_stat, e_rate_ref=1e-3, e_rate=1e-2):
    """
    Return the full emission Probability directly in Log Space.
    ob_stat: Observed Readcounts [2,l] array of 0/1
    """
    e_mat = extended_genotype_emissions (ref_haps, e_rate_ref)
    e_mat = numpy.log (extended_binom_emission (ob_stat=ob_stat, e_mat=e_mat, e_rate=e_rate))
    assert(numpy.max(e_mat) < 0)  # In LOG Space (Assume Error Model)
    return e_mat


def extended_binom_emission (ob_stat, e_mat, e_rate):
    """
    Gives the emission matrix of observed states
    Return emission matrix [n_loci, 4 + 2*n_ref] of each
    ob_stat: [2, n_loci] Matrix with Nr Ref/Alt Reads in Row0/Row1 (!)
    e_mat: Probabilities of genotypes [n_loci, 4 + 2*n_ref+1, 3]
    """

    # What's the probability of observing a dervided read given hidden genotypes 00 01 11
    p_read = numpy.array([e_rate, 0.5, 1 - e_rate])

    # Calculate the Binomial Likelihoods of RC Data
    rc_tot = numpy.sum (ob_stat, axis=0)
    rc_der = ob_stat[1, :]

    prob_binom = binom.pmf (rc_der[:, None], rc_tot[:, None], p_read[None, :])

    # Sum over each of the 3 possible genotypes
    p_full = numpy.sum(e_mat * prob_binom[:, None, :], axis=2)

    return p_full


def extended_genotype_emissions (ref_haps, e_rate_ref):
    """
    Return Emission Matrix, which describes
    probabilities in Genotypes [n_ref+1, n_loci, 3]
    """
    # we already computed it before, but just do again (could be optimized)
    p = numpy.mean (ref_haps, axis=0)

    n_ref = numpy.shape(ref_haps)[0]
    n_loci = numpy.shape(ref_haps)[1]

    # genotypes are 0 - 00, 1 - 10, 2 - 11
    genotypes = [(0,0),(0,1),(1,1)]
    # we have 4 non-roh states, full short ROH, full long ROH
    p_hgeno = -numpy.ones((n_loci, 4 + 2*n_ref, len(genotypes)))

    # Do the 4 non-roh states
    (idx_to_geno, geno_to_idx) = non_roh_state_map()
    assert (len(idx_to_geno) == 4)
    # they all get a one for the right genotype, and 0 otherwise
    for idx in range(len(idx_to_geno)):
        for gdx in range(len(genotypes)):
            # get the genos
            stateGeno = idx_to_geno[idx]
            emitGeno = genotypes[gdx]
            # and compute the probability of the match/missmatch
            p_hgeno[:, idx, gdx] = (1-e_rate_ref) if (stateGeno[0] == emitGeno[0]) else e_rate_ref
            p_hgeno[:, idx, gdx] *= (1-e_rate_ref) if (stateGeno[1] == emitGeno[1]) else e_rate_ref
            # # double-count middle genotype
            if ((emitGeno[0] + emitGeno[1]) == 1):
                tmp = (1-e_rate_ref) if (stateGeno[0] == emitGeno[1]) else e_rate_ref
                tmp *= (1-e_rate_ref) if (stateGeno[1] == emitGeno[0]) else e_rate_ref
                p_hgeno[:, idx, gdx] += tmp


    # Do the copying states (add some error), also twice
    # bounding indicees
    beginShortRoh = 4
    beginLongRoh = 4+n_ref
    # for convenience, transpose the reference
    t_ref_haps = numpy.transpose(ref_haps)
    # one for short ROH beginShortRoh:beginLongRoh
    p_hgeno[:, beginShortRoh:beginLongRoh, 0] = (t_ref_haps == 0) * numpy.power(1 - e_rate_ref,2) + (t_ref_haps == 1) * numpy.power(e_rate_ref,2)
    p_hgeno[:, beginShortRoh:beginLongRoh, 1] = 2*e_rate_ref*(1-e_rate_ref)
    p_hgeno[:, beginShortRoh:beginLongRoh, 2] = (t_ref_haps == 1) * numpy.power(1 - e_rate_ref,2) + (t_ref_haps == 0) * numpy.power(e_rate_ref,2)
    # one for long ROH beginLongRoh:end
    p_hgeno[:, beginLongRoh:, 0] = (t_ref_haps == 0) * numpy.power(1 - e_rate_ref,2) + (t_ref_haps == 1) * numpy.power(e_rate_ref,2)
    p_hgeno[:, beginLongRoh:, 1] = 2*e_rate_ref*(1-e_rate_ref)
    p_hgeno[:, beginLongRoh:, 2] = (t_ref_haps == 1) * numpy.power(1 - e_rate_ref,2) + (t_ref_haps == 0) * numpy.power(e_rate_ref,2)
        
    # Sanity Check if genotype probabilities sum up to (approx.) 1
    assert(numpy.all(numpy.isclose(numpy.sum(p_hgeno, axis=2), 1, rtol=1e-12, atol=1e-12)))
    assert((numpy.min(p_hgeno) >= 0) & (numpy.max(p_hgeno) <= 1))   # Sanity Check

    return p_hgeno

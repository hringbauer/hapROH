"""
Compute the emission probabilities, depending on the observation:
    P( obs_{locus_i, sample_j} | state_i )
"""

import logging
from typing import NewType

import numpy as np
from scipy import stats

from hapROH.classes.genomicData import DataType, GenomicData

logger = logging.getLogger(__name__)

EmissionProba = NewType("EmissionProba", np.ndarray)
"""
(3, nb_snp, nb_samples) array of emission probabilities P( obs_{locus_i, sample_j} | state_i )
The three states are as follows:
    ROH_REF, ROH_ALT, no_ROH
"""


def get_emi_proba(
    genotype_data: GenomicData,
    allele_freq: np.ndarray,
    error_rate: float,
    error_rate_ref: float = 0,
) -> EmissionProba:
    """
    Build emission probabilities for genotype data

    Parameters:
        genotype_data: GenomicData
        allele_freq: np.ndarray of shape (nb_snp), dtype float,
            containing the freq of the alt_allele
        error_rate: float
            the genotyping error
        error_rate_ref: float
            the genotyping error in the reference pannel
    """
    nb_snp, _nb_samples, _ = genotype_data.data.shape
    if nb_snp != len(allele_freq):
        raise ValueError(
            f"`Genotype data` and `allele_freq` contain different number of SNPs ({nb_snp} vs {len(allele_freq)})"
        )
    match genotype_data.datatype:
        case DataType.AD:
            e_proba = _e_proba_from_read_count(
                genotype_data.data, allele_freq, error_rate
            )
        case DataType.GT:
            e_proba = _e_proba_from_GT_count(
                genotype_data.data.sum(axis=2), allele_freq, error_rate
            )
        case DataType.GT_count:
            e_proba = _e_proba_from_GT_count(
                genotype_data.data.squeeze(axis=2), allele_freq, error_rate
            )
        case DataType.PSEUDOHAP:
            e_proba = _e_proba_from_haploid(
                genotype_data.data.squeeze(axis=2), allele_freq, error_rate
            )
        case _:
            raise NotImplementedError(
                f"Emission probabilities not implemented for DataType {genotype_data.datatype}"
            )

    # missing data: uninformative emission => identical for all states
    e_proba[:, genotype_data.get_missing()] = 1

    # error in the reference panel
    d = error_rate_ref * (e_proba[1] - e_proba[0])
    e_proba[0] += d  # => (1-e) * e_proba[0] + e * e_proba[1]
    e_proba[1] -= d  # => (1-e) * e_proba[1] + e * e_proba[0]
    return EmissionProba(e_proba)


def _e_proba_from_haploid(
    genotype_data: np.ndarray, allele_freq: np.ndarray, error_rate: float
) -> EmissionProba:
    """
    Build emission probabilities for haploid data

    Parameters:
        genotype_data: np.ndarray of shape (nb_snp, nb_samples),
            containing 0 (ref) / 1 (alt) / MISSING_VALUE
        allele_freq: np.ndarray of shape (nb_snp), dtype float,
            containing the freq of the alt_allele
        error_rate: float
            the genotyping error
    """
    logger.debug("Computing emission proba from haploid data")
    nb_snp, nb_samples = genotype_data.shape
    allele_freq = allele_freq[
        :, None
    ]  # broadcast to same shape (nb_snp, 1) as genotype_data

    e_mat = np.empty((3, nb_snp, nb_samples), dtype=float)
    e_mat[0, :] = genotype_data == 0  # ROH with ref -> copying state
    e_mat[1, :] = genotype_data == 1  # ROH with alt -> copying state
    e_mat[2, :] = (genotype_data == 0) * (1 - allele_freq) + (
        genotype_data == 1
    ) * allele_freq  # no ROH -> allele freq

    # add error
    e_mat = (1 - error_rate) * e_mat + error_rate * (1 - e_mat)
    return EmissionProba(e_mat)


def _e_proba_from_GT_count(
    genotype_data: np.ndarray, allele_freq: np.ndarray, error_rate: float
) -> EmissionProba:
    """
    Build emission probabilities for diploid GT count

    Parameters:
        genotype_data: np.ndarray of shape (nb_snp, nb_samples),
            containing 0 / 1 / 2 / MISSING_VALUE the nb of alt alleles
        allele_freq: np.ndarray of shape (nb_snp), dtype float,
            containing the freq of the alt_allele
        error_rate: float
            the genotyping error
    """
    logger.debug("Computing emission proba from GT count")
    nb_snp, nb_samples = genotype_data.shape
    allele_freq = allele_freq[
        :, None
    ]  # broadcast to same shape (nb_snp, 1) as genotype_data
    e_mat = np.empty((3, nb_snp, nb_samples), dtype=float)
    e_mat[0, :] = genotype_data == 0  # ROH with ref -> copying state
    e_mat[1, :] = genotype_data == 2  # ROH with alt -> copying state
    e_mat[2, :] = (
        (genotype_data == 0) * (1 - allele_freq) * (1 - allele_freq)
        + (genotype_data == 1) * 2 * (1 - allele_freq) * allele_freq
        + (genotype_data == 2) * allele_freq * allele_freq
    )  # no ROH -> HW
    # add error
    e_mat = (1 - error_rate) * e_mat + error_rate * (1 - e_mat) / 2
    return EmissionProba(e_mat)


def _e_proba_from_read_count(
    genotype_data: np.ndarray, allele_freq: np.ndarray, error_rate: float
) -> EmissionProba:
    """
    Build emission probabilities for read count data

    Parameters:
        genotype_data: np.ndarray of shape (nb_snp, nb_samples, 2),
            containing the nb of REF/ALT reads
        allele_freq: np.ndarray of shape (nb_snp), dtype float,
            containing the freq of the alt_allele
        error_rate: float
            the genotyping error
    """
    logger.debug("Computing emission proba from read count")
    nb_snp, nb_samples, _ = genotype_data.shape

    # P(read count | genotype) follows binomial sampling
    # proba_obs[snp, sample, genotype], genotype = nb of ALT alleles (0, 1, 2)
    nb_ref = genotype_data[:, :, 0]
    nb_reads = genotype_data.sum(axis=2)
    proba_obs = np.empty((nb_snp, nb_samples, 3), dtype=float)
    proba_obs[:, :, 0] = stats.binom.pmf(nb_ref, nb_reads, 1 - error_rate)  # GT = 0,0
    proba_obs[:, :, 1] = stats.binom.pmf(nb_ref, nb_reads, 0.5)  # GT = 0,1
    proba_obs[:, :, 2] = stats.binom.pmf(nb_ref, nb_reads, error_rate)  # GT = 1,1

    # P(genotype | state)
    # proba_geno[snp, genotype, state], state = ROH_REF, ROH_ALT, no_ROH
    proba_geno = np.zeros((nb_snp, 3, 3), dtype=float)
    proba_geno[:, 0, 0] = 1  # ROH with REF -> GT = 0,0
    proba_geno[:, 2, 1] = 1  # ROH with ALT -> GT = 1,1
    # not ROH => HW
    proba_geno[:, 0, 2] = (1 - allele_freq) ** 2
    proba_geno[:, 1, 2] = 2 * allele_freq * (1 - allele_freq)
    proba_geno[:, 2, 2] = allele_freq**2

    # P(read count | state) = sum_genotype P(read count | genotype) * P(genotype | state)
    e_mat = np.einsum("ljg,lgs->slj", proba_obs, proba_geno)
    return EmissionProba(e_mat)


if __name__ == "__main__":
    nb_snp = 10
    nb_samples = 1
    haplo_data = np.random.binomial(1, 0.2, (nb_snp, nb_samples))

    allele_freq = np.random.rand(nb_snp)

    e_mat = _e_proba_from_haploid(haplo_data, allele_freq, 0.05)
    print(haplo_data)
    print(e_mat)

"""
Compute the emission probabilities, depending on the observation:
    P( obs_{locus_i, sample_j} | state_i )
"""

import logging
from typing import NewType

import numpy as np

from hapROH.classes.genomicData import DataType, GenomicData

logger = logging.getLogger(__name__)

# This is only a static-typing alias -> EmissionProba is nothing more than a np.array
EmissionProba = NewType("EmissionProba", np.ndarray)
"""
(3, nb_snp, nb_samples) array of emission probabilities P( obs_{locus_i, sample_j} | state_i )
The three states are as follows:
    ROH_REF, ROH_ALT, no_ROH
"""


def get_emi_proba(
    genotype_data: GenomicData, allele_freq: np.ndarray, error_rate: float
) -> EmissionProba:
    """
    Build emission probabilities for haploid data

    Parameters:
        genotype_data: GenomicData
        allele_freq: np.ndarray of shape (nb_snp), dtype float,
            containing the freq of the alt_allele
        error_rate: float
            the genotyping error
    """
    match genotype_data.datatype:
        case DataType.AD:
            return _e_proba_from_read_count(genotype_data.data, allele_freq, error_rate)
        case DataType.GT:
            return _e_proba_from_GT_count(
                genotype_data.data.sum(axis=2), allele_freq, error_rate
            )
        case DataType.GT_count:
            return _e_proba_from_GT_count(
                genotype_data.data.squeeze(axis=2), allele_freq, error_rate
            )
        case DataType.PSEUDOHAP:
            return _e_proba_from_haploid(
                genotype_data.data.squeeze(axis=2), allele_freq, error_rate
            )


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
    if len(genotype_data.shape) != 2:
        raise ValueError(
            f"Expected genotype data of shape (nb_snp, nb_samples), got {genotype_data.shape}"
        )
    nb_snp, nb_samples = genotype_data.shape
    if nb_snp != len(allele_freq):
        raise ValueError(
            f"`Genotype data` and `allele_freq` cobtain different number of SNPs ({nb_snp} vs {len(allele_freq)})"
        )
    allele_freq = allele_freq[
        :, None
    ]  # broadcast to same shape (nb_snp, 1) as genotype_data
    e_mat = np.ones((3, nb_snp, nb_samples), dtype=float)  # MISSING_VALUE -> 1
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
    if len(genotype_data.shape) != 3 or genotype_data.shape[2] != 1:
        raise ValueError(
            f"Expected genotype data of shape (nb_snp, nb_samples, 1), got {genotype_data.shape}"
        )
    nb_snp, nb_samples, _ = genotype_data.shape
    if nb_snp != len(allele_freq):
        raise ValueError(
            f"`Genotype data` and `allele_freq` cobtain different number of SNPs ({nb_snp} vs {len(allele_freq)})"
        )
    allele_freq = allele_freq[
        :, None
    ]  # broadcast to same shape (nb_snp, 1) as genotype_data
    e_mat = np.ones((3, nb_snp, nb_samples), dtype=float)  # MISSING_VALUE -> 1
    e_mat[0, :] = genotype_data == 0  # ROH with ref -> copying state
    e_mat[1, :] = genotype_data == 2  # ROH with alt -> copying state
    e_mat[2, :] = (
        (genotype_data == 0) * (1 - allele_freq) * (1 - allele_freq)
        + (genotype_data == 1) * 2 * (1 - allele_freq) * allele_freq
        + (genotype_data == 2) * allele_freq * allele_freq
    )  # no ROH -> HW
    # add error
    e_mat = (1 - error_rate) * e_mat + error_rate * (1 - e_mat)
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
    raise NotImplementedError("Emission for read count not implemented yet")


if __name__ == "__main__":
    nb_snp = 10
    nb_samples = 1
    haplo_data = np.random.binomial(1, 0.2, (nb_snp, nb_samples))

    allele_freq = np.random.rand(nb_snp)

    e_mat = _e_proba_from_haploid(haplo_data, allele_freq, 0.05)
    print(haplo_data)
    print(e_mat)

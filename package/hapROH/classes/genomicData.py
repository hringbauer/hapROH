from __future__ import annotations

import logging
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import StrEnum

import h5py
import numpy as np
import numpy.typing as npt
import pandas as pd

logger = logging.getLogger(__name__)

MISSING_VAL = -1

################################
# Class describing genomic data
################################


class DataType(StrEnum):
    AD = "readcount"  # two columns with read_counts
    GT = "GT"  # two columns with 0 (ref) / 1 (alt) / MISSING_VALUE
    PSEUDOHAP = "haploid"  # one column with 0 (ref) / 1 (alt) / MISSING_VALUE
    GT_count = (
        "diploid_gt"  # one column with 0 / 1 / 2 / MISSING_VALUE the nb of alt alleles
    )


@dataclass
class GenomicData:
    data: np.ndarray  # shape (nb_snp, nb_samples, 1|2), dtype int or float
    datatype: DataType

    def flip_data(self, idx_flipped: npt.NDArray[np.bool_]) -> None:
        """Flip the ref and alt allele at the given positions (in place modification)."""
        match self.datatype:
            case DataType.AD | DataType.GT:
                self.data[idx_flipped] = self.data[idx_flipped][:, :, ::-1]
            case DataType.GT_count:
                mask_missing = self.data[idx_flipped] == MISSING_VAL
                self.data[idx_flipped] = np.where(
                    mask_missing, MISSING_VAL, 2 - self.data[idx_flipped]
                )
            case DataType.PSEUDOHAP:
                mask_missing = self.data[idx_flipped] == MISSING_VAL
                self.data[idx_flipped] = np.where(
                    mask_missing, MISSING_VAL, 1 - self.data[idx_flipped]
                )
            case _:
                raise NotImplementedError(
                    f"flip_data not implemented for DataType {self.datatype}"
                )

    def to_GT_count(
        self, allele_freq: None | np.ndarray = None, error_rate: float = 1e-3
    ) -> GenomicData:
        """Transform the data into (diploid) genotype counts.

        For DataType.AD, genotypes are called with a simple Bayesian model
        combining a Hardy-Weinberg prior derived from `allele_freq`
        and a per read sequencing `error_rate`

        Parameters
        ----------
        allele_freq : np.ndarray
            Alt allele frequency at each SNP, shape (nb_snp,).
            Only used when self.datatype is DataType.AD.
        error_rate : float
            Sequencing error rate for each read.
            Only used when self.datatype is DataType.AD.
        """
        match self.datatype:
            case DataType.AD:
                logger.warning(
                    "Calling genotypes from readcounts using bayesian model. This is not recommended for low coverage data. In that case, prefer e_model=haploid"
                )
                if allele_freq is None:
                    raise ValueError(
                        "Please provide an allele frequency to call genotype from allele depth"
                    )
                n_ref = self.data[:, :, 0].astype(np.float64)
                n_alt = self.data[:, :, 1].astype(np.float64)
                depth = n_ref + n_alt
                is_missing = depth == 0

                # prior based on HW-equilibrium using allele_freq
                p = np.clip(allele_freq, 1e-6, 1 - 1e-6)  # shape (nb_snp,)
                log_prior = np.stack(
                    [
                        2 * np.log(1 - p),
                        np.log(2) + np.log(p) + np.log(1 - p),
                        2 * np.log(p),
                    ],
                    axis=-1,
                )  # (nb_snp, 3)

                # likelihood under a binomial model to observe (n_ref, n_alt)
                alt_frac = (
                    error_rate,
                    0.5,
                    1 - error_rate,
                )  # expected alt-read fraction for GT = 0, 1, 2
                log_lik = np.stack(
                    [n_alt * np.log(f) + n_ref * np.log(1 - f) for f in alt_frac],
                    axis=-1,
                )  # (nb_snp, nb_samples, 3)

                log_post = log_lik + log_prior[:, np.newaxis, :]
                new_data = np.argmax(log_post, axis=-1).astype(np.int8)[..., np.newaxis]
                new_data[is_missing] = MISSING_VAL

            case DataType.GT:
                mask_missing = np.any(self.data == MISSING_VAL, axis=2)
                new_data = np.sum(self.data, axis=2, keepdims=True).astype(np.int8)
                new_data[mask_missing] = MISSING_VAL

            case DataType.GT_count:
                new_data = self.data.copy()

            case DataType.PSEUDOHAP:
                logger.warning("Trying to convert haploid to diploid data. ")
                mask_missing = self.data == MISSING_VAL
                new_data = 2 * self.data
                new_data[mask_missing] = MISSING_VAL
                # raise ValueError("Cannot recover diploid genotype counts from pseudo-haploid data ")
            case _:
                raise NotImplementedError(
                    f"to_GT_count not implemented for DataType {self.datatype}"
                )

        return GenomicData(new_data, DataType.GT_count)

    def to_pseudo_haploid(self, seed: None | int = None) -> GenomicData:
        """Transform the data into pseudo-haploid data."""
        rng = np.random.default_rng(seed)

        match self.datatype:
            case DataType.AD:
                depth = np.sum(self.data, axis=2)
                is_missing = depth == 0
                alt_prob = self.data[:, :, 1] / np.maximum(
                    1, depth
                )  # avoid division by 0
                new_data = rng.binomial(1, alt_prob).astype(np.int8)[..., np.newaxis]
                new_data[is_missing] = MISSING_VAL

            case DataType.GT:
                choice_idx = rng.integers(0, 2, size=self.data.shape[:2])
                new_data = np.take_along_axis(
                    self.data, choice_idx[..., np.newaxis], axis=2
                ).astype(np.int8)

            case DataType.GT_count:
                missing = self.data == MISSING_VAL
                het_pos = self.data == 1
                new_data = np.where(self.data == 2, 1, 0).astype(np.int8)
                random_het = rng.integers(0, 2, size=self.data.shape).astype(np.int8)
                new_data[het_pos] = random_het[het_pos]
                new_data[missing] = MISSING_VAL

            case DataType.PSEUDOHAP:
                new_data = self.data.copy()
            case _:
                raise NotImplementedError(
                    f"to_pseudo_haploid not implemented for DataType {self.datatype}"
                )

        return GenomicData(new_data, DataType.PSEUDOHAP)

    def downsample(
        self, target_depth: float = 1, seed: None | int = None
    ) -> GenomicData:
        """Downsample the data to a given depth.
        Works only if datatype is AD"""
        if not self.datatype == DataType.AD:
            raise ValueError(f"Cannot downsample datatype {self.datatype}, only AD")
        mean_depth_per_sample = np.mean(np.sum(self.data, axis=2), axis=0)
        if np.any(mean_depth_per_sample <= target_depth):
            low_mask = mean_depth_per_sample <= target_depth
            low_samples = [
                f"{i} ({d:.3f}x)"
                for i, d in zip(low_mask.nonzero(), mean_depth_per_sample[low_mask])
            ]
            raise ValueError(
                f"Target depth {target_depth:.3} is higher than actual mean depth for samples {low_samples}. Cannot downsample !"
            )
        p = target_depth / mean_depth_per_sample
        new_data = np.random.default_rng(seed).binomial(
            self.data, p[np.newaxis, :, np.newaxis]
        )
        return GenomicData(new_data, DataType.AD)


class GenomicDataFile(ABC):
    @abstractmethod
    def get_iids(self) -> npt.NDArray[np.bytes_]:
        """Return a string array with the sample names present in the file"""

    @abstractmethod
    def get_snp(self) -> pd.DataFrame:
        """
        Return a pandas dataframe with the following columns:
            pos (int), map (float, in Morgans), ref (U1), alt (U1), chrom (int, optional)
        """

    @abstractmethod
    def get_data(
        self,
        idx_snp: None | npt.NDArray[np.bool_] = None,
        idx_iid: None | npt.NDArray[np.bool_] = None,
    ) -> GenomicData:
        """
        Return the genomic data
        Args:
            idx_snp: if provided, boolean mask indicating which positions to keep
            idx_iid: if provided, boolean mask indicating which individuals to keep
        """

    def get_idx_iids(self, iids: str | list[str]) -> npt.NDArray[np.bool_]:
        """Return a boolean mask corresponding to the position of given individual(s) in the dataset"""
        if isinstance(iids, str):
            iids = [iids]
        idx_iids = np.isin(self.get_iids(), iids)
        if sum(idx_iids) != len(iids):
            raise ValueError(
                f"Found {sum(idx_iids)} matching individuals instead of the {len(iids)} queried ones"
            )
        return idx_iids

    @classmethod
    def load_genetic_file(cls, path: str) -> GenomicDataFile:
        """
        Load genomic data from a given file.
        File type is detected automatically (currently supported: eigenstrat, hdf5)
        """
        ext_hdf5 = [".hdf5", ".h5"]
        ext_eigenstrat = [".snp", ".geno", ".ind"]
        if os.path.isfile(path):
            path_prefix, ext = os.path.splitext(path)
            if ext in ext_hdf5:
                return Hdf5File(path)
            elif ext in ext_eigenstrat:
                return EigenstratFile(path_prefix)
            else:
                raise ValueError(
                    f"File extension {ext} not recognised. File should be of type hdf5 (.hdf5|.h5) or eigenstrat"
                )
        for ext in ext_hdf5:
            if os.path.isfile(path + ext):
                return Hdf5File(path + ext)
        for ext in ext_eigenstrat:
            if os.path.isfile(path + ext):
                return EigenstratFile(path)
        raise FileNotFoundError(f"File {path} not found")


################################
# Eigenstrat implementation
################################


def load_geno_unpacked(
    path: str,
    idx_snp: None | npt.NDArray[np.bool_] = None,
    idx_iid: None | npt.NDArray[np.bool_] = None,
) -> np.ndarray:
    """Load genotype from an unpacked eigenstrat geno file."""
    columns = None if idx_iid is None else idx_iid.nonzero()[0].tolist()
    geno = np.genfromtxt(path, usecols=columns, dtype="i1", delimiter=1)
    geno = 2 - geno  # replace the nb of ref alleles by the nb of alt alleles
    geno[geno == 2 - 9] = MISSING_VAL
    row_idx = slice(None) if idx_snp is None else idx_snp
    return geno[row_idx, :, np.newaxis]  # shape (nb_snp, nb_samples, 1)


def load_geno_packed(
    path: str,
    idx_snp: None | npt.NDArray[np.bool_] = None,
    idx_iid: None | npt.NDArray[np.bool_] = None,
) -> np.ndarray:
    """Load genotype from a packed eigenstrat geno file."""
    # read header to get number of samples and of snp
    with open(path, "rb") as f:
        header_bytes = f.read(48)
    header = header_bytes.decode("ascii", errors="ignore").split()
    magic, nb_samples_str, nb_snp_str = header[0], header[1], header[2]
    assert magic == "GENO", f"Unexpected magic string {magic!r}, expected 'GENO'"
    nb_samples = int(nb_samples_str)
    nb_snp = int(nb_snp_str)

    rlen = max(48, (2 * nb_samples + 7) // 8)  # bytes per record

    raw = np.fromfile(path, dtype=np.uint8)
    expected_size = rlen * (nb_snp + 1)
    assert raw.size == expected_size, (
        f"file size {raw.size} != expected {expected_size} based on the number of samples and SNPs"
    )
    raw = raw[rlen:].reshape(nb_snp, rlen)

    bits = np.unpackbits(raw, axis=1)
    bits = bits[:, : 2 * nb_samples].reshape(nb_snp, nb_samples, 2)
    geno = (2 * bits[:, :, 0] + bits[:, :, 1]).astype(np.int8)  # in 0..3

    geno = 2 - geno  # replace the nb of ref alleles by the nb of alt alleles
    geno[geno == 2 - 3] = MISSING_VAL
    if idx_snp is not None:
        geno = geno[idx_snp]
    if idx_iid is not None:
        geno = geno[:, idx_iid]
    return geno[..., np.newaxis]


class EigenstratFile(GenomicDataFile):
    path_prefix: str

    def __init__(self, path_prefix: str):
        self.path_prefix = path_prefix

    def get_iids(self) -> npt.NDArray[np.bytes_]:
        return np.loadtxt(self.path_prefix + ".ind", usecols=0, dtype="U50")

    def get_snp(self) -> pd.DataFrame:
        return pd.read_csv(
            self.path_prefix + ".snp",
            header=None,
            sep=r"\s+",
            names=["SNP", "chrom", "map", "pos", "ref", "alt"],
        )

    def get_data(
        self,
        idx_snp: None | npt.NDArray[np.bool_] = None,
        idx_iid: None | npt.NDArray[np.bool_] = None,
    ) -> GenomicData:
        logger.info(f"Loading data from {self.path_prefix}")
        geno_path = self.path_prefix + ".geno"
        with open(geno_path, "rb") as f:
            magic = f.read(4)

        if magic == b"GENO":
            data = load_geno_packed(geno_path, idx_snp, idx_iid)
        else:
            data = load_geno_unpacked(geno_path, idx_snp, idx_iid)

        datatype = DataType.GT_count
        # if np.all(data[:,:,0]!=1):            # search for heterozygotes
        #     datatype = DataType.PSEUDOHAP
        #     data[data == 2] = 1

        logger.info(f"Loaded {data.shape} SNPs, of dtype {data.dtype}")
        return GenomicData(data, datatype)


################################
# Hdf5 implementation
################################


class Hdf5File(GenomicDataFile):
    path: str
    filter_biallelic_snp: bool
    mask_snp: None | npt.NDArray[np.bool_] = (
        None  # changed only if filter_biallelic_snp is True
    )

    def __init__(self, path: str, filter_biallelic_snp: bool = True):
        self.path = path
        self.filter_biallelic_snp = filter_biallelic_snp

    def get_iids(self) -> npt.NDArray[np.bytes_]:
        with h5py.File(self.path, "r") as h5_file:
            samples = np.array(h5_file["samples"]).astype("U50")
        return samples

    def get_snp(self) -> pd.DataFrame:
        with h5py.File(self.path, "r") as h5_file:
            variants = h5_file["variants"]
            pos = np.array(variants["POS"]).astype(int)
            map = np.array(variants["MAP"]).astype(float)
            ref = np.array(variants["REF"]).astype("U1")
            alt = np.array(variants["ALT"]).astype("U1")
            chrom = (
                np.array(variants["CHROM"]).astype(int) if "CHROM" in variants else None
            )

        df_snp = pd.DataFrame({"pos": pos, "map": map, "ref": ref, "alt": alt})
        if chrom is not None:
            df_snp["chrom"] = chrom

        if self.filter_biallelic_snp:
            # filter only biallelic variants and reshape alt into a 1d-array if needed
            mask_biallelic = np.ones(len(alt), dtype=bool)
            if alt.ndim == 2:
                mask_biallelic = np.all(alt[:, 1:] == "", axis=1)
                alt = alt[:, 0]

            # filter only SNP
            bases = np.array(["A", "T", "G", "C"])
            mask_snp = np.isin(ref, bases) & np.isin(alt, bases)

            mask = mask_biallelic & mask_snp
            logger.info(f"Kept {np.sum(mask)}/{len(mask)} biallelic SNP sites")
            df_snp = df_snp[mask]
            self.mask_snp = mask
        return df_snp

    def get_data(
        self,
        idx_snp: None | npt.NDArray[np.bool_] = None,
        idx_iid: None | npt.NDArray[np.bool_] = None,
    ) -> GenomicData:
        logger.info(f"Loading data from {self.path}")
        if self.filter_biallelic_snp:
            if self.mask_snp is None:
                self.get_snp()
            assert self.mask_snp is not None, (
                "self.mask_snp is set by self.get_snp(), when self.filter_bialleleic_snp is True "
            )
            row_idx = self.mask_snp
            if idx_snp is not None:
                row_idx = idx_snp & self.mask_snp
        else:
            row_idx = slice(None) if idx_snp is None else idx_snp

        column_idx = slice(None) if idx_iid is None else idx_iid

        with h5py.File(self.path, "r") as h5_file:
            calldata = h5_file["calldata"]
            if "AD" in calldata.keys():
                data = calldata["AD"][:][row_idx][:, column_idx]
                datatype = DataType.AD
                if "GT" in calldata.keys():
                    logger.warning(
                        f"{self.path} contains both fields AD and GT. Only AD is loaded, GT is ignored"
                    )
            elif "GT" in calldata.keys():
                data = calldata["GT"][:][row_idx][:, column_idx]
                datatype = DataType.GT
            else:
                raise ValueError(
                    f"Found neither AD nor GT field in the hdf5 file {self.path}"
                )
        if datatype == DataType.GT:
            # check if pseudohaploid and reshape if necessary
            if len(data.shape) == 2:
                data = data[..., np.newaxis]
            if len(data.shape) == 3:
                if data.shape[2] == 1:
                    datatype = DataType.GT_count
                    # if np.all(data != 2):                   # search for heterozygotes
                    #     datatype = DataType.PSEUDOHAP
                elif data.shape[2] == 2:
                    datatype = DataType.GT
                    # if np.all(data[:,:,0]==data[:,:,1]):    # search for heterozygotes
                    #     datatype = DataType.PSEUDOHAP
                    #     data = data[:,:,:1]
                else:
                    raise ValueError(
                        f"Expected data of shape (nb_snp, nb_samples, 1|2), not {data.shape}"
                    )
            else:
                raise ValueError(
                    f"Expected data of shape (nb_snp, nb_samples, 1|2), not {data.shape}"
                )
        logger.info(f"Loaded {data.shape} SNPs, of dtype {data.dtype}")
        return GenomicData(data, datatype)


########################################################
# Utility methods to work with genomic data
########################################################


def get_snp_intersection(
    snp_sample: pd.DataFrame, snp_ref: pd.DataFrame, chrom: None | int = None
) -> tuple[npt.NDArray[np.bool_], npt.NDArray[np.bool_], npt.NDArray[np.bool_]]:
    """Get the indices of the intersecting positions in the two SNP sets.

    Returns (sample_idx, ref_idx, mismatching_sample_idx) as numpy arrays.
        - `sample_idx`/`ref_idx` are boolean masks with intersecting positions from `snp_sample`/`snp_ref` where REF/ALT are identical or flipped
        - `mismatching_sample_idx` is a boolean mask with positions from snp_sample[sample_idx] where REF/ALT are flipped compared to snp_ref[ref_idx]
    """
    # restrict to a given chromosome, if necessary
    idx_chr_sample = np.ones(len(snp_sample), dtype=bool)
    idx_chr_ref = np.ones(len(snp_ref), dtype=bool)
    if chrom is not None:
        if "chrom" in snp_sample.columns:
            idx_chr_sample = snp_sample["chrom"] == chrom
            if sum(idx_chr_sample) == 0:
                raise ValueError(
                    f"Chromosome {chrom} not found in data (chromosomes present: {sorted(snp_sample['chrom'].unique())})"
                )
        else:
            logger.debug("Field `chrom` not found in `snp_sample`")
        if "chrom" in snp_ref.columns:
            idx_chr_ref = snp_ref["chrom"] == chrom
            if sum(idx_chr_ref) == 0:
                raise ValueError(
                    f"Chromosome {chrom} not found in data (chromosomes present: {sorted(snp_ref['chrom'].unique())})"
                )
        else:
            logger.debug("Field `chrom` not found in `snp_ref`")

    # keep track of the original row order, to know where each position came from
    new_snp_sample = snp_sample.assign(_sample_idx=np.arange(len(snp_sample)))
    new_snp_ref = snp_ref.assign(_ref_idx=np.arange(len(snp_ref)))
    # mask other chromosomes
    new_snp_sample = new_snp_sample[idx_chr_sample]
    new_snp_ref = new_snp_ref[idx_chr_ref]
    merged = new_snp_ref.merge(
        new_snp_sample[["pos", "_sample_idx", "ref", "alt"]],
        on="pos",
        how="inner",
    ).sort_values("pos")

    # check if ref and alt are flipped
    matching_SNP = (merged["ref_x"] == merged["ref_y"]) & (
        merged["alt_x"] == merged["alt_y"]
    )  # boolean array, len(merged)
    flipped_SNP = (merged["ref_x"] == merged["alt_y"]) & (
        merged["alt_x"] == merged["ref_y"]
    )  # boolean array, len(merged)
    mismatching_SNP = ~matching_SNP & ~flipped_SNP  # boolean array, len(merged)

    logger.info(
        f"{len(merged)}/{len(snp_ref)} SNP found in intersection, of which {flipped_SNP.sum()} flipped REF/ALT and {mismatching_SNP.sum()} mismatching REF/ALT"
    )

    sample_mask = np.zeros(len(snp_sample), dtype=bool)
    ref_mask = np.zeros(len(snp_ref), dtype=bool)
    sample_mask[merged["_sample_idx"][~mismatching_SNP]] = (
        True  # boolean array, len(snp_sample)
    )
    ref_mask[merged["_ref_idx"][~mismatching_SNP]] = True  # boolean array, len(snp_ref

    return sample_mask, ref_mask, flipped_SNP[~mismatching_SNP].to_numpy()


def get_rmap(
    df_snp: pd.DataFrame, min_gap: float = 1e-10, max_gap: float = np.inf
) -> np.ndarray:
    """Return the genetic distance [in Morgan] between locci, clipping values to the desired interval"""
    gen_pos = df_snp["map"]
    assert gen_pos.is_monotonic_increasing, (
        "SNP positions must be sorted in ascending order"
    )
    if gen_pos.max() > 20:
        logger.debug("Converting from centimorgans to morgans")
        gen_pos /= 100
    r_map = gen_pos[1:].to_numpy() - gen_pos[:-1].to_numpy()
    logger.info(f"Minimum Genetic Map: {gen_pos.min()} Morgan")
    logger.info(f"Maximum Genetic Map: {gen_pos.max()} Morgan")
    logger.info(f"Gaps bigger than 0.1 cM: {(r_map > 0.001).sum()}")
    logger.info(f"Maximum Gap: {r_map.max() * 100:.4f} cM")
    logger.info(f"Clipping gaps to range: {100 * min_gap:.3f} - {100 * max_gap:.3f} cM")
    return np.clip(r_map, min_gap, max_gap)


if __name__ == "__main__":
    nb_snp = 10
    nb_iids = 3
    read_counts = np.arange(2 * nb_snp * nb_iids).reshape((nb_snp, nb_iids, 2))
    gt = np.random.binomial(1, 0.5, 2 * nb_snp * nb_iids).reshape((nb_snp, nb_iids, 2))

    ad_data = GenomicData(read_counts, DataType.AD)
    gt_data = GenomicData(gt, DataType.GT)

    ### Check that conversion return expected shape
    # AD -> downsampled AD
    ad_data_down = ad_data.downsample(2)
    assert ad_data_down.data.shape == (nb_snp, nb_iids, 2), (
        f"{ad_data_down.data.shape} != {(nb_snp, nb_iids, 2)}"
    )

    # AD -> GT_count
    allele_freq = np.full(nb_snp, 0.5)
    gtc_data = ad_data.to_GT_count(allele_freq)
    assert gtc_data.data.shape == (nb_snp, nb_iids, 1), (
        f"{gtc_data.data.shape} != {(nb_snp, nb_iids, 1)}"
    )

    # AD -> pseudo-haploid
    hap_data = ad_data.to_pseudo_haploid()
    assert hap_data.data.shape == (nb_snp, nb_iids, 1), (
        f"{hap_data.data.shape} != {(nb_snp, nb_iids, 1)}"
    )

    # GT -> GT_count
    gtc_data = gt_data.to_GT_count(allele_freq)
    assert gtc_data.data.shape == (nb_snp, nb_iids, 1), (
        f"{gtc_data.data.shape} != {(nb_snp, nb_iids, 1)}"
    )

    # GT -> pseudo-hap
    hap_data = gt_data.to_pseudo_haploid()
    assert hap_data.data.shape == (nb_snp, nb_iids, 1), (
        f"{hap_data.data.shape} != {(nb_snp, nb_iids, 1)}"
    )

    # GT_count -> pseudo-hap
    hap_data = gtc_data.to_pseudo_haploid()
    assert hap_data.data.shape == (nb_snp, nb_iids, 1), (
        f"{hap_data.data.shape} != {(nb_snp, nb_iids, 1)}"
    )

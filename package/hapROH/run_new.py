import logging
import os
from typing import Literal

import numpy as np

from hapROH.classes.genomicData import (
    DataType,
    GenomicDataFile,
    get_rmap,
    get_snp_intersection,
)
from hapROH.classes.HMM import HMM
from hapROH.classes.postprocessing import postproces

logger = logging.getLogger(__name__)


def callROH_chr(
    path_sample: str,
    path_ref: str,
    chrom: int,
    iids: None | str | list[str] = None,
    folder_out: str = "",
    r_in: float = 1,
    r_out: float = 20,
    r_jump: float = 300,
    error_rate: float = 0.01,
    e_model: Literal["readcounts", "diploid_gt", "haploid"] = "haploid",
    downsampling: None | float = None,
    cutoff_post: float = 0.999,
    snps_extend: int = 0,
    max_gap: float = 0.005,
    min_len1: float = 0.04,
    min_len2: float = 0.02,
    min_len_final: float = 0.04,
    logfile: None | str = None,
    loglevel: int = 0,
    backend: Literal["python", "numba", "cython"] = "cython",
):
    """Call runs of homozygosity (ROH) for one chromosome using an HMM.

    Loads sample and reference genotype data for a single chromosome, intersects
    their SNP positions (flipping REF/ALT alleles where needed), builds a genetic
    map for the intersected SNPs, optionally downsamples and/or transforms the
    sample data according to the chosen emission model, runs the ROH HMM to
    compute posterior state probabilities, and writes the resulting ROH calls
    and auxiliary per-SNP data to disk for each requested individual.

    Args:
        path_sample: Path to the sample genotype file (Eigenstrat or HDF5) to call ROH on.
        path_ref: Path to the reference panel genotype file used for allele
            frequencies/haplotypes.
        chrom: Chromosome number to process.
        iids: Individual ID(s) to call ROH for. If None, all individuals in the
            sample file are used. A single string is treated as one IID.
        folder_out: Base output directory. Results for each individual are written
            to ``folder_out/<iid>/chr<chrom>/``.
        r_in: HMM transition rate into the ROH (homozygous) state.
        r_out: HMM transition rate out of the ROH state.
        r_jump: HMM jump rate between two distinct ROH states.
        error_rate: Genotyping/sequencing error rate used by the emission model.
        e_model: Model to use for computing the emission probabilities. One of:
            - "readcounts": use raw read-count (AD) data; requires the sample
              data to contain an 'AD' field.
            - "diploid_gt": use diploid genotype counts computed against the
              reference allele frequencies; not valid for haploid input data.
            - "haploid": pseudo-haploid calls.
        downsampling: If provided, depth to downsample the sample data
            to before calling ROH. Only valable if input data is AD
        cutoff_post: Cutoff used when calling the ROH segments from the posterior probability
        snps_extend: Number of SNPs added to elongate ROH blocks (before merging).
        max_gap: Maximum gap (in Morgans) between two adjacent ROH for them to be merged.
        min_len1 and min_len2: Minimum lengths (in Morgans) required for two adjacent ROH to be merged.
        min_len_final: Minimum length (in Morgans) for segments to appear in the final dataset.
        logfile: Path to a file to write log output to. If None, logs go to the default stream handler.
        loglevel: Verbosity level for the hapROH logger (0=WARNING, 1=INFO, 2 or higher=DEBUG).
        backend: Which backend to the for the forward-backward algorithm.

    Returns:
        np.ndarray: The posterior probability array returned by the HMM, with
        shape (nb_state=nb_ref_haplo+1, n_snps, n_individuals).

    Written files:
        - logfile: if logfile is provided
        For each individual in `iids`, writes the following files under
        ``folder_out/<iid>/chr<chrom>/``:
            - roh.csv: called ROH segments (with columns for iid and chrom).
            - posterior0.csv: per-SNP posterior probability of the ROH state.
            - pos.csv: physical positions of the intersected SNPs.
            - map.csv: genetic map positions of the intersected SNPs.
            - readcounts.csv / gt_count.csv: per-SNP genotype/read-count data,
              depending on the sample data's datatype (not written for
              pseudo-haploid data).
    """
    level = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}.get(
        loglevel, logging.DEBUG
    )
    logging.basicConfig(
        format="%(asctime)s.%(msecs)03d [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        filename=logfile,
        level=logging.WARNING,
    )
    logging.captureWarnings(True)
    logging.getLogger("hapROH").setLevel(level)

    logger.info(f"Starting callROH_chr on chromosome {chrom} and iids {iids}")
    logger.info(f"Sample file: {path_sample}")
    logger.info(f"Reference file: {path_ref}")
    logger.info(f"HMM parameters: r_in={r_in}, r_out={r_out}, r_jump={r_jump}")
    logger.info(f"Emission model={e_model}, error_rate={error_rate}")
    logger.info(
        f"Merging parameters: max_gap={max_gap}, (min_len1, min_len2)={(min_len1, min_len2)}"
    )
    logger.info(f"Filtering parameters: min_len_final={min_len_final}")

    if min_len_final < min(min_len1, min_len2):
        logger.warning(
            f"min_len_final < min(min_len1,min_len2). By design, all segments shorter than {min(min_len1, min_len2)} will still be removed."
        )

    ### Preload the files
    file_sample = GenomicDataFile.load_genetic_file(path_sample)
    file_ref = GenomicDataFile.load_genetic_file(path_ref)

    ### Filter individuals
    if iids is None:
        iids = file_sample.get_iids().astype(str).tolist()
    elif isinstance(iids, str):
        iids = [iids]
    assert isinstance(iids, list)
    idx_iids = file_sample.get_idx_iids(iids)
    # TODO: if wanted, filter iids in file_ref

    ### Load the SNPs and compute intersection
    logger.info("Loading SNP and computing intersection")
    snp_sample = file_sample.get_snp()
    snp_ref = file_ref.get_snp()
    idx_snp_sample, idx_snp_ref, idx_flipped_sample = get_snp_intersection(
        snp_sample, snp_ref, chrom
    )
    logger.info(
        f"Found {sum(idx_snp_sample)} intersecting SNP, of which {sum(idx_flipped_sample)} flipped REF/ALT"
    )

    ### Get the recombinaton map + restict SNP to intersecting positions
    df_snp = snp_ref.iloc[idx_snp_ref]
    r_map = get_rmap(df_snp)

    ### Load the genomic data at the intersecting positions
    logger.info("Loading genotype data")
    data_sample = file_sample.get_data(idx_snp_sample, idx_iids)
    data_ref = file_ref.get_data(idx_snp_ref)
    logger.info("Done loading genotype data")

    data_sample.flip_data(idx_flipped_sample)

    ### Preprocess the data
    data_sampl_pp = data_sample
    if downsampling is not None:
        data_sampl_pp = data_sampl_pp.downsample(downsampling)
    match e_model:
        case "readcount":
            if data_sampl_pp.datatype != "readcount":
                raise ValueError(
                    "Cannot use e_model='readcount' if imput data does not contain a field 'calldata/AD'"
                )
        case "diploid_gt":
            if data_sampl_pp.datatype == "haploid":
                raise ValueError(
                    "Cannot use e_model='diploid_gt' if imput data is haploid."
                )
            data_sampl_pp = data_sampl_pp.to_GT_count(
                allele_freq=data_ref.data.mean(axis=(1, 2))
            )
        case "haploid":
            data_sampl_pp = data_sampl_pp.to_pseudo_haploid()
        case _:
            raise ValueError(
                f"Invalid option e_model='{e_model}'. Valid options are None | 'readcounts' | 'diploid_gt' | 'haploid'"
            )

    ### Initialise the HMM
    logger.info("Initialising HMM")
    hmm = HMM(data_sampl_pp, data_ref, r_map, r_in, r_out, r_jump, error_rate)
    logger.info("Done initialising HMM")

    ### Compute the posterior probability
    logger.info("Computing posterior probabilities")
    post_pb = hmm.calc_posterior_proba(backend)
    logger.info("Done computing posterior probabilities")

    ### Postprocess and save the results
    logger.info("Saving result")
    for idx, iid in enumerate(iids):
        folder_out_iid = os.path.join(folder_out, iid, "chr" + str(chrom), "")
        logger.info(f"Writing individual {iid} to {folder_out_iid}")
        if not os.path.isdir(folder_out_iid):
            os.makedirs(folder_out_iid)

        # Create and save ROH dataframe
        df_roh_iid = postproces(
            post_pb[0, :, idx],
            df_snp,
            cutoff_post,
            snps_extend,
            max_gap,
            min_len1,
            min_len2,
        )
        df_roh_iid = df_roh_iid[df_roh_iid["lengthM"] >= min_len_final]
        df_roh_iid["iid"] = iid
        df_roh_iid["ch"] = chrom
        df_roh_iid.to_csv(folder_out_iid + "roh.csv", index=False)

        # Save SNP info
        np.savetxt(
            folder_out_iid + "posterior0.csv",
            post_pb[0, :, idx],
            delimiter=",",
            fmt="%f",
        )
        np.savetxt(folder_out_iid + "pos.csv", df_snp["pos"], delimiter=",", fmt="%f")
        np.savetxt(folder_out_iid + "map.csv", df_snp["map"], delimiter=",", fmt="%f")

        # Save GT info along results for latter plotting
        match data_sample.datatype:
            case DataType.AD:
                hap = np.zeros((2, len(df_snp["map"])), dtype="uint8")  # dummy data
                np.savetxt(
                    folder_out_iid + "readcounts.csv",
                    data_sample.data[:, idx],
                    delimiter=",",
                    fmt="%f",
                )
            case DataType.GT:
                hap = data_sample.data[:, idx].T
            case DataType.GT_count:
                hap = np.empty((2, len(df_snp["map"])), dtype="uint8")
                hap[0] = np.where(data_sample.data[:, idx].squeeze(axis=-1) == 2, 1, 0)
                hap[1] = np.where(data_sample.data[:, idx].squeeze(axis=-1) != 0, 1, 0)
            case DataType.PSEUDOHAP:
                hap = np.tile(data_sample.data[:, idx], (2, 1))
            case _:
                raise NotImplementedError("Case not implemented")
        np.savetxt(folder_out_iid + "hap.csv", hap, delimiter=",", fmt="%f")

    return post_pb

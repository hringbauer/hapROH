import logging
from typing import List

import numpy as np
import pandas as pd

from .classes.genomicData import DataType, GenomicDataFile, get_snp_intersection
from .classes.HMM import HMM
from .utils.miscellanious import print_memory_usage

logger = logging.getLogger(__name__)

def get_rmap(df_snp:pd.DataFrame, min_gap:float=1e-10, max_gap:float=np.inf) -> np.ndarray:
    """Return the genetic distance [in Morgan] between locci, clipping values to the desired interval"""
    gen_pos = df_snp["map"]
    assert gen_pos.is_monotonic_increasing, "SNP positions must be sorted in ascending order"
    if gen_pos.max() > 20:
        logger.debug(f"Converting from centimorgans to morgans")
        gen_pos /= 100
    r_map = gen_pos[1:].to_numpy() - gen_pos[:-1].to_numpy()
    logger.info(f"Minimum Genetic Map: {gen_pos.min()} Morgan")
    logger.info(f"Maximum Genetic Map: {gen_pos.max()} Morgan")
    logger.info(f"Gaps bigger than 0.1 cM: {(r_map > 0.001).sum()}")
    logger.info(f"Maximum Gap: {r_map.max() * 100:.4f} cM")
    logger.info(f"Clipping gaps to range: {100*min_gap:.3f} - {100*max_gap:.3f} cM")
    return np.clip(r_map, min_gap, max_gap)

def callROH_chr(path_sample:str, path_ref:str, chrom:int, iids:None|str|List[str]=None,
                r_in:float=1, r_out:float=20, r_jump: float=300, error_rate:float=0.01,
                e_model:None|str="haploid", downsampling:None|float=None,
                logfile:None|str=None, loglevel:int=0
                ):
    level = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}.get(loglevel, logging.DEBUG)
    logging.basicConfig(
        format="%(asctime)s.%(msecs)03d [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        filename=logfile,
        level=level
    )
    logging.captureWarnings(True)

    logger.info(f"Starting callROH_chr on chromosome {chrom}")
    print_memory_usage(logger)

    ### Preload the files
    file_sample = GenomicDataFile.load_genetic_file(path_sample)
    file_ref = GenomicDataFile.load_genetic_file(path_ref)

    ### Filter individuals
    idx_iids = None if iids is None else file_sample.get_idx_iids(iids)
    # TODO: if wanted, filter iids in file_ref

    ### Load the SNPs and compute intersection
    logger.info(f"Loading SNP and computing intersection")
    snp_sample = file_sample.get_snp()
    snp_ref = file_ref.get_snp()
    idx_snp_sample, idx_snp_ref, idx_flipped_sample = get_snp_intersection(snp_sample, snp_ref, chrom)
    logger.info(f"Found {sum(idx_snp_sample)} intersecting SNP, of which {sum(idx_flipped_sample)} flipped REF/ALT")
    print_memory_usage(logger)

    ### Get the recombinaton map
    df_snp = snp_ref.iloc[idx_snp_ref]
    r_map = get_rmap(df_snp)
    print_memory_usage(logger)

    ### Load the genomic data at the intersecting positions
    logger.info(f"Loading genotype data")
    data_sample = file_sample.get_data(idx_snp_sample, idx_iids)
    print_memory_usage(logger)
    data_ref = file_ref.get_data(idx_snp_ref)
    logger.info(f"Done loading genotype data")
    print_memory_usage(logger)

    data_sample.flip_data(idx_flipped_sample)
    print_memory_usage(logger)

    ### If wanted: preprocess the data
    if downsampling is not None:
        data_sample.downsample(downsampling)
    if e_model is not None:
        match e_model:
            case "readcount":
                if data_sample.datatype != "readcount":
                    raise ValueError(f"Cannot use readcount model if provided data does not contain a field calldata/AD")
            case "diploid_gt":
                data_sample.to_GT_count(allele_freq = data_ref.data.mean(axis=(1,2)))
            case "haploid":
                data_sample.to_pseudo_haploid()
    print_memory_usage(logger)

    if not data_ref.datatype == DataType.GT:
        raise ValueError(f"Reference pannel should contain GT but contains {data_ref.datatype}")

    ref_panel = data_ref.data.reshape(data_ref.data.shape[0], -1)   # (nb_snp, nb_samples, 2) -> (nb_snp, 2*nb_samples)
    print_memory_usage(logger)

    ### Load the actual HMM
    logger.info("Initialising HMM")
    hmm = HMM(data_sample, ref_panel, r_map,
                    r_in, r_out, r_jump, error_rate)
    logger.debug("Done initialising HMM")

    # logger.info("Computing posterior probabilities")
    # post_pb = hmm.calc_posterior_proba()
    # logger.info("Done computing posterior probabilities")

    return None
    return post_pb
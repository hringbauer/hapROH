import logging, os
from typing import List

import numpy as np

from hapROH.classes.genomicData import DataType, GenomicDataFile, get_snp_intersection, get_rmap
from hapROH.classes.HMM import HMM
from hapROH.classes.postprocessing import postproces
from hapROH.utils.miscellanious import print_memory_usage

logger = logging.getLogger(__name__)

def callROH_chr(path_sample:str, path_ref:str, chrom:int, iids:None|str|List[str]=None,
                folder_out:str="",
                r_in:float=1, r_out:float=20, r_jump: float=300, error_rate:float=0.01,
                e_model:None|str="haploid", downsampling:None|float=None,
                cutoff_post:float=0.999,
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

    logger.info(f"Starting callROH_chr on chromosome {chrom} and iids {iids}")
    logger.info(f"Sample file: {path_sample}")
    logger.info(f"Reference file: {path_ref}")

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
    logger.info(f"Loading SNP and computing intersection")
    snp_sample = file_sample.get_snp()
    snp_ref = file_ref.get_snp()
    idx_snp_sample, idx_snp_ref, idx_flipped_sample = get_snp_intersection(snp_sample, snp_ref, chrom)
    logger.info(f"Found {sum(idx_snp_sample)} intersecting SNP, of which {sum(idx_flipped_sample)} flipped REF/ALT")

    ### Get the recombinaton map + restict SNP to intersecting positions
    df_snp = snp_ref.iloc[idx_snp_ref]
    r_map = get_rmap(df_snp)

    ### Load the genomic data at the intersecting positions
    logger.info(f"Loading genotype data")
    data_sample = file_sample.get_data(idx_snp_sample, idx_iids)
    data_ref = file_ref.get_data(idx_snp_ref)
    logger.info(f"Done loading genotype data")

    data_sample.flip_data(idx_flipped_sample)

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

    ### Initialise the HMM
    logger.info("Initialising HMM")
    hmm = HMM(data_sample, data_ref, r_map,
                    r_in, r_out, r_jump, error_rate)
    logger.debug("Done initialising HMM")

    ### Compute the posterior probability
    logger.debug("Computing posterior probabilities")
    post_pb = hmm.calc_posterior_proba()
    logger.debug("Done computing posterior probabilities")

    ### Postprocess and save the results
    logger.debug(f"Saving result")
    for idx, iid in enumerate(iids):
        folder_out_iid = os.path.join(folder_out, iid, "chr" + str(chrom), "")
        logger.debug(f"Writing individual {iid} to {folder_out_iid}")
        if not os.path.isdir(folder_out_iid):
            os.makedirs(folder_out_iid)

        # Creating ROH dataframe
        df_roh_iid = postproces(post_pb[0,:,idx], df_snp, cutoff_post)
        df_roh_iid["iid"] = iid
        df_roh_iid["chrom"] = chrom
        df_roh_iid.to_csv(folder_out_iid+"roh.csv", index=False)

        # Saving SNP info
        np.savetxt(folder_out_iid+"posterior0.csv", post_pb[0, :, idx], delimiter=",",  fmt='%f')
        np.savetxt(folder_out_iid+"pos.csv", df_snp["pos"], delimiter=",",  fmt='%f')
        np.savetxt(folder_out_iid+"map.csv", df_snp["map"], delimiter=",",  fmt='%f')

    return post_pb
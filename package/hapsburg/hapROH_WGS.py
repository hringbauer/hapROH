import argparse
import h5py
import pandas as pd
import os
from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind_lowmem
from hapsburg.PackagesSupport.loadEigenstrat.saveHDF5 import mpileup2hdf5
from hapsburg.PackagesSupport.parallel_runs.helper_functions import multi_run
from multiprocessing import set_start_method

import numpy as np

def main():
    set_start_method("spawn")
    parser = argparse.ArgumentParser(description='Running hapROH on WGS BAM file.')
    parser.add_argument('-m', action="store", dest="m", type=str, required=True,
                        help="path to folder containing the pileup file. The pileup file should be named as {args.iid}.chr*.mpileup")
    parser.add_argument('-r', action="store", dest="ref", type=str, required=True,
                        help="path to reference panel")
    parser.add_argument('--meta', action="store", dest="meta", type=str, required=True, 
                        default="/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv",
                        help="path to metadata of the reference panel")
    parser.add_argument('-i', action="store", dest="iid", type=str, required=False, default="",
                        help="IID of the target individual. If unspecified, will use the prefix of the bam file.")
    parser.add_argument('-q', action="store", dest="q", type=int, required=False, default=25, 
                        help="Minimum mapping quality for reads to be considered.")
    parser.add_argument('-Q', action="store", dest="Q", type=int, required=False, default=30, 
                        help="Minimum base quality.")
    parser.add_argument('-p', action="store", dest="processes", type=int, required=False, default=1,
                        help="Number of processes to use.")
    parser.add_argument('--out', action="store", dest="out", type=str, required=False, default="",
                        help="folder to store output files")
    args = parser.parse_args()

    # check if the reference panel has binary-encoded genotypes
    # with h5py.File(args.ref, 'r') as f:
    #     if not ('calldata' in f.keys() and 'GTbinary' in f['calldata'].keys()):
    #         raise RuntimeError(f'The reference panel {args.r} does not have binary-encoded genotypes. Please use a reference panel with the field calldata/binaryGT.')
    
    iid = args.iid
    prms = [ [os.path.join(args.m, f'{iid}.chr{ch}.mpileup'),
        args.ref + str(ch) + ".hdf5",
        iid, -np.inf, np.inf, os.path.join(args.out, 'hdf5'), False] for ch in range(1, 23)]
    results = multi_run(mpileup2hdf5, prms, args.processes)
    
    err = np.mean(np.array([err for err, _, _, _ in results]))
    numSitesCovered = sum([n for _, n, _, _ in results])
    totNumSites = sum([n for _, _, n, _ in results])
    print(f'Estimated genotyping error: {err}')
    frac = numSitesCovered/totNumSites
    print(f'Fraction of sites covered by >=1 read: {frac}')
    
    if frac >= 0.7:
        downsample = 1.0
    else:
        downsample = False

    ################################### call ROH ##################################
    ## first, run without contamination
    df = hapsb_ind_lowmem(iid, chs=range(1,23), 
        path_targets_prefix = os.path.join(args.out, 'hdf5'),
        h5_path1000g = args.ref, meta_path_ref = args.meta,
        folder_out=os.path.join(args.out, 'hapROH'), prefix_out="",
        e_model="readcount", p_model="HDF5", post_model="Standard",
        processes=args.processes, delete=False, output=True, save=True, save_fp=False, 
        n_ref=2504, diploid_ref=True, exclude_pops=[], readcounts=True, random_allele=False, downsample=downsample, 
        c=0.0, roh_min_l_final=0.04, roh_in=1, roh_out=20, roh_jump=300, e_rate=err, e_rate_ref=1e-3, 
        logfile=True, combine=True, file_result="_roh_full.csv")
    
if __name__ == "__main__":
    main()
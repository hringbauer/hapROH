"""
Run Eigenstrat inference on the cluster.
Called with array jobs from sbatch
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import numpy as np
import os as os
import sys as sys
import pandas as pd
import socket

#### 1) Set the Path to default HAPSBURG
path = "/project2/jnovembre/hringbauer/HAPSBURG/"  # The Path on Midway Cluster
os.chdir(path)
from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind  # Need this import

#########################################################

def get_iid_path(i, reps=20, down_sampling_covs = np.linspace(0.2, 1.0, 9),
                 base_path="./Data/SA_1240kHDF5/Ust_Ishim_rep/downsample_ph_r"):
    """Get the Individual IID"""
    batch = int(np.floor(i/reps))
    rep = i%reps
    
    path_hd = base_path + str(rep) + ".h5"
    c = down_sampling_covs[batch]
    iid = f"{c:.4f}_r{rep}" 
    return iid, path_hd

#########################################################
#########################################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Script needs argument (indiviual i)")
    run_nr = int(sys.argv[1]) # The Parameter passed to the Python Script from outside

    iid, path_target = get_iid_path(i=run_nr)
    
    hapsb_ind(iid=iid, chs=range(1, 23), 
              path_targets=path_target, # The path before the .ind, .snp, .geno
              h5_path1000g='./Data/1000Genomes/HDF5/1240kHDF5/all1240int8/chr', 
              meta_path_ref='./Data/1000Genomes/Individuals/meta_df_all.csv', 
              folder_out="./Empirical/1240k/SA_Readcounts/Ust_Ishim_rep/", prefix_out='', 
              e_model='haploid', p_model='MosaicHDF5', 
              post_model='Standard', processes=1, delete=False, output=True, save=True, 
              save_fp=False, n_ref=2504, exclude_pops=[], readcounts=False, random_allele=True, 
              roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.0, 
              cutoff_post=0.999, max_gap=0, roh_min_l=0.01, 
              logfile=True, combine=True, file_result='_roh_full.csv')
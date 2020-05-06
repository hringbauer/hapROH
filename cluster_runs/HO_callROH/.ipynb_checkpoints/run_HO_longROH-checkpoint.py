"""
Run Human Origin inference for long ROH on cluster.
Called with array jobs from sbatch
@ Author: Harald Ringbauer, 2020, All rights reserved
"""

import numpy as np
import os as os
import sys as sys
import pandas as pd

#### 1) Set the Path to default HAPSBURG
path = "/project2/jnovembre/hringbauer/HAPSBURG/"  # The Path on Midway Cluster
os.chdir(path)

### 1.5) Do some follow up imports
sys.path.append("./package/")  
from hapsburg.PackagesSupport.hapsburg_run import hapsb_chrom, hapsb_ind

#########################################################

def get_iid_from_df(i, path_df="./Empirical/HO/combined_roh05.csv", 
                    sep="\t", cm=12):
    """Return ith iid to run"""
    df = pd.read_csv(path_df, sep=sep)
    idx = df[f"sum_roh>{cm}"]>0
    print(f"Indivdiuals with ROH>{cm}: {np.sum(idx)}/{len(idx)}")
    iids = df["iid"][idx].values
    return iids[i]


#########################################################
#########################################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Script needs argument (indiviual i)")
    
    run_nr = int(sys.argv[1]) # The Parameter passed to the Python Script from outside
    iid=get_iid_from_df(run_nr)
    
    print(f"Runnding haploid mode for {iid}...")
    hapsb_ind(iid=iid, chs=range(1, 23), 
              path_targets= "./Data/Marcus2019_1240k/sardinia_hapsburg.h5", 
              h5_path1000g= './Data/1000Genomes/HDF5/1240kHDF5/all1240int8/chr', 
              meta_path_ref= './Data/1000Genomes/Individuals/meta_df_all.csv', 
              folder_out = './Empirical/HO/long_hap/', 
              prefix_out = '', e_model="haploid", p_model="SardHDF5", post_model='Standard', 
              processes=1, delete=False, output=True, save=True, save_fp=False, 
              n_ref=2504, exclude_pops=[], readcounts=False, 
              random_allele=True, roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.0, 
              cutoff_post=0.999, max_gap=0, roh_min_l=0.01, logfile=True, 
              combine=True, file_result='_roh_full.csv')
    
    print(f"Runnding diploid mode for {iid}...")
    hapsb_ind(iid=iid, chs=range(1, 23), 
              path_targets= "./Data/Marcus2019_1240k/sardinia_hapsburg.h5", 
              h5_path1000g='./Data/1000Genomes/HDF5/1240kHDF5/all1240int8/chr', 
              meta_path_ref='./Data/1000Genomes/Individuals/meta_df_all.csv', 
              folder_out='./Empirical/HO/long_dip/', 
              prefix_out='', e_model="diploid_gt", p_model="SardHDF5", post_model='Standard', 
              processes=1, delete=False, output=True, save=True, save_fp=False, n_ref=2504, exclude_pops=[], readcounts=False, 
              random_allele=False, roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.0, 
              cutoff_post=0.999, max_gap=0, roh_min_l=0.01, logfile=True, combine=True, file_result='_roh_full.csv')
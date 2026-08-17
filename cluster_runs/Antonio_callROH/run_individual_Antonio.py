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

### 1.5) Do follow up imports
sys.path.append("./PackagesSupport/")
from hapsburg_run import hapsb_chrom, hapsb_ind

### 2) Load the Meta Data
#meta_path = "./Data/Antonio/meta_processed.csv"  # All Individuals to run
meta_path = "./PackagesSupport/cluster_runs/Antonio_callROH/rerun.csv"  # For rerun
meta_df = pd.read_csv(meta_path, sep="\t")
print(f"Loaded {len(meta_df)} Ancients")


#########################################################

def get_iid_from_i(df, i):
    """Get the Individual IID #i from Dataframe df_anno when filtering for good coverage"""    
    if i<0 or i>=len(df):    # Sanity Check
        raise RuntimeError(f"Index {i} out of Range of High Coverage ancients.") 
    iid = df["iid"].values[i]
    print(f"Running individual {iid} \n")
    return iid


#########################################################
#########################################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Script needs single argument (indiviual i)")
        
    run_nr = int(sys.argv[1]) # The Parameter passed to the Python Script from outside
    iid = get_iid_from_i(meta_df, run_nr)
    
    hapsb_ind(iid=iid, chs=range(1,23), processes=1, 
          h5_path_targets = "./Data/Antonio/rmpr_unigeno_1240k.hdf5",
          base_out_folder="./Empirical/1240k/Antonio/",
          e_model="readcount", p_model="MosaicHDF5", n_ref=2504,
          delete=False, logfile=True, combine=True)
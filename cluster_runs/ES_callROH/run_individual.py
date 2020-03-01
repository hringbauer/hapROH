"""
Run Eigenstrat inference on the cluster.
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


####################################################################
####################################################################
### Functions


def load_eigenstrat_anno(path="./Data/ReichLabEigenstrat/Raw/v37.2.1240K.clean4.anno", 
                         anc_only=True, min_snps=400000):
    """Load annotated Eigenstrat (from D. Reich's group).
    anc_only: Return only the ancients with age>0"""
    df_anno = pd.read_csv(path, sep="\t", engine="python")
    coverage = pd.to_numeric(df_anno["Coverage"], errors='coerce')
    df_anno["coverage"]=coverage

    # Convert the Ages as well
    ages = df_anno["Average of 95.4% date range in calBP (defined as 1950 CE)  "]
    df_anno["ages"] = pd.to_numeric(ages, errors='coerce')  #

    ### Convert Longitude and Latitude
    lat = df_anno["Lat."]
    lon = df_anno["Long."]
    df_anno["lat"] = pd.to_numeric(lat, errors='coerce')
    df_anno["lon"] = pd.to_numeric(lon, errors='coerce')
    
    df_anc = df_anno[df_anno["ages"]>0]

    print(f"Loaded {len(df_anc)} / {len(df_anno)} ancient Indivdiuals Anno File.")
    
    if anc_only:
        df_anno=df_anc
        
    df_anno = df_anno[df_anno["SNPs hit on autosomes"]>min_snps]
    print(f"Loaded {len(df_anno)} Individuals with >{min_snps} SNPs covered")
    return df_anno

def get_iid_from_df(df, i, id_col="Instance ID"):
    """Get the Individual IID"""
    if i<0 or i>=len(df):    # Sanity Check
        raise RuntimeError(f"Index {i} out of Range of High Coverage ancients.") 
    iid = df[id_col].values[i]
    return iid

#########################################################
#########################################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Script needs argument (indiviual i)")
    
    run_nr = int(sys.argv[1]) # The Parameter passed to the Python Script from outside
    df_anno = load_eigenstrat_anno(min_snps=400000)
    iid = get_iid_from_df(df_anno, run_nr, id_col="Instance ID")
    
    hapsb_ind(iid, chs=range(1, 23), processes=1, delete=False, output=True, 
               save=True, save_fp=False, n_ref=2504, exclude_pops=[], 
               e_model='haploid', p_model='EigenstratPacked', readcounts=False, 
               destroy_phase=True, post_model='Standard', 
               path_targets='./Data/ReichLabEigenstrat/Raw/v37.2.1240K', 
               h5_path1000g='./Data/1000Genomes/HDF5/1240kHDF5/all1240int8/chr', 
               meta_path_ref='./Data/1000Genomes/Individuals/meta_df_all.csv', 
               base_out_folder='./Empirical/Eigenstrat/Reichall/final/', prefix_out='', 
               roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.0, max_gap=0, 
               cutoff=0.999, l_cutoff=0.01, logfile=True, combine=True, file_name='_roh_full.csv')
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

### 1.5) Do some follow up imports
sys.path.append("./Python3/")  
sys.path.append("./PackagesSupport/parallel_runs/")
from hmm_inference import HMM_Analyze   # The HMM core object
from helper_functions import prepare_path, combine_individual_data

### 2) Load the Data
meta_path = "./Data/Marcus2019_1240k/meta_rev_final.csv"
anc_sardind= 85
anc_ind =  1087
meta_df = pd.read_csv(meta_path)
df_anno = meta_df[anc_sardind:anc_ind]
print(f"Loaded {len(df_anno)} Ancients")


####################################################################
####################################################################
### Functions

def analyze_chromosome_ph(iid, ch=3, n_ref=503, save=True, save_fp=False, exclude_pops=[], 
                          base_out_folder="./Empirical/Eigenstrat/Reichall/", prefix_out="",
                          roh_in=100, roh_out=100, roh_jump=300, e_rate=0.01, e_rate_ref=0.01, 
                          max_gap=0, logfile=True):
    """Run the analysis for one individual and chromosome on Pseudohaploid data
    Wrapper for HMM Class. Takes 13 Parameters"""
    
    ### The folder on what to run the Data on (PERMANENTLY set here to fixed loaction)
    ## What target Files to run on:
    h5_path_targets = "./Data/Marcus2019_1240k/mod_reich_sardinia_ancients_rev_mrg_dedup_3trm_anno.h5"
    meta_path_targets = "./Data/Marcus2019_1240k/meta_rev_final.csv"
    
    ## Reference Files:
    h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/all1240/chr" 
    meta_path_ref = "./Data/1000Genomes/Individuals/meta_df_all.csv"
    
    ### Create Folder if needed, and pipe output if wanted
    path_out = prepare_path(base_out_folder, iid, ch, prefix_out, logfile=logfile)
    
    hmm = HMM_Analyze(cython=2, p_model="SardHDF5", e_model="haploid", post_model="Standard",
                      manual_load=True, save=save, save_fp=save_fp)
    

    ### Load and prepare the pre-processing Model
    hmm.load_preprocessing_model()              # Load the preprocessing Model
    hmm.p_obj.set_params(readcounts = False, destroy_phase=True, base_out_folder=base_out_folder, 
                         prefix_out_data=prefix_out, excluded=exclude_pops)   
    
    ### Set paths of reference and target to run on. DELETE ref when run for with European Reference!
    hmm.p_obj.set_params(h5_path1000g = h5_path1000g, meta_path_ref = meta_path_ref,
                         h5_path_targets = h5_path_targets, meta_path_targets=meta_path_targets)
    
    hmm.load_data(iid=iid, ch=ch, n_ref=n_ref)  # Load the actual Data
    hmm.load_secondary_objects()

    ### Set the Parameters
    hmm.e_obj.set_params(e_rate = e_rate, e_rate_ref = e_rate_ref)
    hmm.t_obj.set_params(roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump)
    hmm.post_obj.set_params(max_gap=max_gap)
    
    #hmm.calc_viterbi_path(save=save)           # Calculate the Viterbi Path.
    hmm.calc_posterior(save=save)              # Calculate the Posterior.
    hmm.post_processing(save=save)             # Do the Post-Processing.
                         
#########################################################
def analyze_individual_ph(iid, chs=range(1,23), n_ref=2504, save=True, save_fp=False, 
                          exclude_pops=[], base_out_folder="./Empirical/Eigenstrat/Reichall/", 
                          prefix_out="", roh_in=100, roh_out=100, roh_jump=300, e_rate=0.001, 
                          e_rate_ref=0.001, max_gap=0, logfile=True, output=True, processes=5, delete=True):
    """Analyze a full single individual in a parallelized fasion. Run all Chromosome analyses in parallel
    Wrapper for analyze_chromosome_gt.
    logfile: Whether to use a logfile
    output: Whether to print general Output"""
                            
    if output == True:
        print(f"Doing Individual {iid}...")
    
    ### Prepare the Parameters for that Indivdiual
    prms = [[iid, ch, n_ref, save, save_fp, exclude_pops, base_out_folder, prefix_out,
         roh_in, roh_out, roh_jump, e_rate, e_rate_ref, max_gap, logfile] for ch in chs] 
                            
    ### Run the analysis in parallel
        # Run the analysis for all Parameters
    for prm in prms:
        assert(len(prm)==15)  # Sanity Check
        analyze_chromosome_ph(*prm)
                            
    ### Merge results for that Individual
    combine_individual_data(base_out_folder, iid=iid, delete=delete, chs=chs)                  
    return

#########################################################

def get_iid_from_i(df_anno, i, min_cov):
    """Get the Individual IID #i from Dataframe df_anno when filtering for good coverage"""
    df_t = df_anno[(df_anno["mean_cov"] > min_cov) & (df_anno["include_alt"] > 0)]
    print(f"Nr. high coverage ancients: {len(df_t)}")
    
    if i<0 or i>=len(df_t):    # Sanity Check
        raise RuntimeError(f"Index {i} out of Range of High Coverage ancients.") 
    iid = df_t["iid"].values[i]
    print(f"Picking individual {iid} \n")
    return iid


#########################################################
#########################################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Script needs single argument (indiviual i)")
        
    run_nr = int(sys.argv[1]) # The Parameter passed to the Python Script from outside
    iid = get_iid_from_i(df_anno, run_nr, min_cov=0.5)
    analyze_individual_ph(iid=iid, chs=range(1, 23), n_ref=2504, save=True, save_fp=False, exclude_pops=[],
                          base_out_folder="./Empirical/1240k/MarcusAncsPH/", prefix_out="",
                          roh_in=100, roh_out=100, roh_jump=300, e_rate=0.01,
                          e_rate_ref=0.01, max_gap=0.0, logfile=False, output=True, delete=True)

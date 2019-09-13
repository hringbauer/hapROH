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


def load_eigenstrat_anno(path="./Data/ReichLabEigenstrat/Raw/v37.2.1240K.clean4.anno", anc_only=True):
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
    return df_anno

#df_anno = load_eigenstrat_anno()
df_anno = pd.read_csv("./PackagesSupport/cluster_runs/ES_callROH/rerun.csv")  # For the rerun of Shotgun Individuals



####################################################################
### Functions

def analyze_chromosome_es(iid, ch=3, n_ref=503, save=True, save_fp=False, exclude_pops=[], 
                          base_out_folder="./Empirical/Eigenstrat/Reichall/", prefix_out="",
                          roh_in=100, roh_out=100, roh_jump=300, e_rate=0.01, e_rate_ref=0.001, 
                          max_gap=0, logfile=True):
    """Run the analysis for one individual and chromosome on eigenstrat data
    Wrapper for HMM Class. Takes 13 Parameters"""
    
    ### The folder on what to run the Data on (PERMANENTLY set here to fixed loaction)
    ## What Eigenstrat File to run on:
    es_target_path="./Data/ReichLabEigenstrat/Raw/v37.2.1240K"
    
    ## Reference Files:
    h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/all1240/chr" 
    meta_path_ref = "./Data/1000Genomes/Individuals/meta_df_all.csv"
    
    ### Create Folder if needed, and pipe output if wanted
    path_out = prepare_path(base_out_folder, iid, ch, prefix_out, logfile=logfile)
    
    hmm = HMM_Analyze(cython=2, p_model="Eigenstrat", e_model="haploid", post_model="Standard",
                      manual_load=True, save=save, save_fp=save_fp)

    ### Load and prepare the pre-processing Model
    hmm.load_preprocessing_model()              # Load the preprocessing Model
    hmm.p_obj.set_params(es_target_path=es_target_path, readcounts = False, destroy_phase=False,
                base_out_folder=base_out_folder, prefix_out_data=prefix_out, excluded=exclude_pops)   
    
    ### Set to run with full 1000G reference. DELETE when run for with European Reference!!
    hmm.p_obj.set_params(h5_path1000g = h5_path1000g, meta_path_ref = meta_path_ref)
    
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
def analyze_individual_es(iid, chs=range(1,23), n_ref=2504, save=True, save_fp=False, 
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
        analyze_chromosome_es(*prm)
                            
    ### Merge results for that Individual
    combine_individual_data(base_out_folder, iid=iid, delete=delete, chs=chs)                  
    return

#########################################################

def get_iid_from_i(df_anno, i, min_cov=0.5):
    """Get the Individual IID"""
    df_t = df_anno[df_anno["coverage"] > min_cov] # Extract high coverage individuals
    if i<0 or i>=len(df_t):    # Sanity Check
        raise RuntimeError(f"Index {i} out of Range of High Coverage ancients.") 
    iid = df_t["Instance ID"].values[i]
    return iid


#########################################################
#########################################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Script needs argument (indiviual i)")
    run_nr = int(sys.argv[1]) # The Parameter passed to the Python Script from outside

    iid = get_iid_from_i(df_anno, run_nr, min_cov=0.5)
    analyze_individual_es(iid=iid, chs=range(1, 23), n_ref=2504, save=True, save_fp=False, exclude_pops=[],
                          base_out_folder="./Empirical/Eigenstrat/Reichall/", prefix_out="",
                          roh_in=100, roh_out=100, roh_jump=300, e_rate=0.01,
                          e_rate_ref=0.001, max_gap=0.0, logfile=False, output=True, delete=True)

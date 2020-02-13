"""
Helper Function to run HAPSBURG on a full individual and full reference Dataset.
Have function for running on single chromsome, with ALL keywords as variables.
@Author: Harald Ringbauer, September 2019
All rights reserved
"""

import numpy as np
import os as os
import sys as sys
import multiprocessing as mp
import pandas as pd
import socket

# Assume that now we are in the root directory
#sys.path.append("./XX/") Go to base directory of HAPSBURG
from hapsburg.hmm_inference import HMM_Analyze   # The HMM core object

#sys.path.append("./PackagesSupport/parallel_runs/")
#from helper_functions import prepare_path, multi_run, combine_individual_data
from hapsburg.PackagesSupport.parallel_runs.helper_functions import prepare_path, multi_run, combine_individual_data


def hapsb_chrom(iid, ch=3, save=True, save_fp=False, n_ref=2504, exclude_pops=[],
                e_model="readcount", p_model="MosaicHDF5", readcounts=True, destroy_phase=True,
                post_model="Standard", path_targets = "./Data/SA_1240kHDF5/IPK12.h5",
                h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/all1240/chr", 
                meta_path_ref = "./Data/1000Genomes/Individuals/meta_df_all.csv",
                base_out_folder="./Empirical/Eigenstrat/Reichall/test/", prefix_out="",
                roh_in=100, roh_out=100, roh_jump=300, e_rate=0.01, e_rate_ref=0.01, 
                max_gap=0, logfile=True):
    """Run Hapsburg analysis for one chromosome on eigenstrat data
    Wrapper for HMM Class. Takes 20 Parameters"""
    
    ### Create Folder if needed, and pipe output if wanted
    path_out = prepare_path(base_out_folder, iid, ch, prefix_out, logfile=logfile)
    
    hmm = HMM_Analyze(cython=2, p_model=p_model, e_model=e_model,
                      manual_load=True, save=save, save_fp=save_fp)

    ### Load and prepare the pre-processing Model
    hmm.load_preprocessing_model()              # Load the preprocessing Model
    hmm.p_obj.set_params(readcounts = readcounts, destroy_phase=destroy_phase,
                base_out_folder=base_out_folder, prefix_out_data=prefix_out, excluded=exclude_pops)
    
    ### Set the paths to ref & target
    hmm.p_obj.set_params(h5_path1000g = h5_path1000g, path_targets = path_targets, 
                         meta_path_ref = meta_path_ref, n_ref=n_ref)
    hmm.load_data(iid=iid, ch=ch)  # Load the actual Data
    hmm.load_secondary_objects()

    ### Set the Parameters
    hmm.e_obj.set_params(e_rate = e_rate, e_rate_ref = e_rate_ref)
    hmm.t_obj.set_params(roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump)
    hmm.post_obj.set_params(max_gap=max_gap)
    
    ### hmm.calc_viterbi_path(save=save)           # Calculate the Viterbi Path.
    hmm.calc_posterior(save=save)              # Calculate the Posterior.
    hmm.post_processing(save=save)             # Do the Post-Processing.
         

#########################################################
### Run Hapsburg for one Individual (wrap for Chr.)

def hapsb_ind(iid, chs=range(1,23), processes=1, delete=False, output=True, save=True, save_fp=False, n_ref=2504, 
              exclude_pops=[], e_model="readcount", p_model="MosaicHDF5", readcounts=True, destroy_phase=False,
              post_model="Standard", path_targets = "./Data/SA_1240kHDF5/IPK12.h5",
              h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/all1240/chr", 
              meta_path_ref = "./Data/1000Genomes/Individuals/meta_df_all.csv",
              base_out_folder="./Empirical/Eigenstrat/Reichall/test/", prefix_out="",
              roh_in=100, roh_out=100, roh_jump=300, e_rate=0.01, e_rate_ref=0.01, 
              max_gap=0, logfile=True, combine=True, file_name="_roh_full.csv"):
    """Analyze a full single individual in a parallelized fasion. Run all Chromosome analyses in parallel
    Wrapper for hapsb_chrom
    logfile: Whether to use a logfile
    output: Whether to print general output.
    delete: Whether to delete the output-folder after successfully combining output file
    combine: Whether to combine output into one .csv per Individual (with suffix file_name)
    DEFAULT is with parameters to do Readcounts from a Reference HDF5 (and 1000G Ref)"""
                            
    if output:
        print(f"Doing Individual {iid}...")
    
    ### Prepare the Parameters for that Indivdiual
    prms = [[iid, ch, save, save_fp, n_ref, exclude_pops, e_model, p_model, readcounts, destroy_phase,
            post_model, path_targets, h5_path1000g, meta_path_ref, base_out_folder, prefix_out,
            roh_in, roh_out, roh_jump, e_rate, e_rate_ref, max_gap, logfile] for ch in chs]
    assert(len(prms[0])==23)   # Sanity Check
                            
    ### Run the analysis in parallel
    multi_run(hapsb_chrom, prms, processes = processes)
                            
    ### Merge results for that Individual
    if combine:
        if output:
            print(f"Combining Information for {len(chs)} Chromosomes...")
        combine_individual_data(base_out_folder, iid=iid, delete=delete, chs=chs, prefix_out=prefix_out, file_name=file_name)
    if output:
        print(f"Run finished successfully!")
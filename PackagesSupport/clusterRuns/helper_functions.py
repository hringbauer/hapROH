"""
Helper Functions for Notebook Runs on Cluster
@ Author: Harald Ringbauer, 2019, All rights reserved
"""
import numpy as np
import pandas as pd
import sys
import os

def prepare_path(path_mosaic, iid, ch, prefix_out, logfile=True):
    """Prepare the path and pipe printing for one Individual.
    Create Path if not already existing.
    logfile: Whether to pipe output to log-file"""       
    path_out = os.path.join(path_mosaic, "output", iid, "chr" + str(ch), prefix_out, "")
    if not os.path.exists(path_out):
            os.makedirs(path_out)
    ### Active LOG FILE if needed
    if logfile == True:
        #path_log = path_log + "hmm_run_log.txt"
        path_log = os.path.join(path_out, "hmm_run_log.txt")
        print(f"Set Output Log path: {path_log}")
        sys.stdout = open(path_log, 'w') 
    return path_out


def multi_run(fun, prms, processes = 4):
    """Implementation of running in Parallel.
    fun: Function
    prms: The Parameter Files
    processes: How many Processes to use"""
    print(f"Running {len(prms)} jobs in parallel.")
    
    with mp.Pool(processes = processes) as pool:
        results = pool.starmap(fun, prms)

def split_up_roh_df(base_path, path_out, iid):
    """Splits up the ROH-dataframe.
    base_path: Where to find roh_info.csv
    path_out: Where to save roh_gt to
    iid: Which Individual to extract from roh_info.csv"""
    #path = base_path + "roh_info.csv"
    path = os.path.join(base_path, "roh_info.csv")
    dft = pd.read_csv(path, sep="\t")  # Load the Meta File

    save_df = dft[dft["iid"] == iid]
    save_path = os.path.join(path_out, "roh_gt.csv")
    save_df.to_csv(save_path, sep="\t", index=False)
    #print(f"Saved to {save_path}")
    return

def combine_individual_data(base_path, iid, delete=False, chs=range(1,23), prefix_out=""):
    """Function to merge data from one Individual Analysis (all Chromosome)
    chs: Which Chromosomes to combine"
    delete: Whether to delete individual folder and contents after combining."""
    
    full_df_vec =[]  # The full dataframe of inferred ROH blocks
    
    ### Walk through Chromosomes and combine the Dataframes
    for ch in chs:
        path_roh = os.path.join(base_path, str(iid), "chr"+str(ch), prefix_out, "roh.csv") 
        df_temp = pd.read_csv(path_roh, sep=",")
        full_df_vec.append(df_temp)
        
    full_df = pd.concat(full_df_vec)
        
    ### Save to Path:
    path_save = os.path.join(base_path, str(iid) + "_roh_full.csv")
    full_df.to_csv(path_save, index=False)
    
    ### Delete files in folder if need
    if delete == True:
        for ch in chs:
            path_folder = os.path.join(base_path, str(iid), "chr"+str(ch), prefix_out, "") 
            
            for root, _, files in os.walk(path_folder):
                for file in files:
                    os.remove(os.path.join(root, file))
            os.rmdir(path_folder) # Remove the Chromosome Folders
        os.rmdir(os.path.join(base_path, str(iid), ""))  # Remove the Individual Folder
    
    return full_df







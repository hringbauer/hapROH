"""
Helper Functions for Notebook Runs on Cluster
@ Author: Harald Ringbauer, 2019
"""

import numpy as np
import pandas as pd
import sys
import os
import multiprocessing as mp

def prepare_path(base_path, iid, ch, prefix_out, logfile=True):
    """Prepare the path and pipe printing for one Individual.
    Create Path if not already existing.
    logfile: Whether to pipe output to log-file"""       
    path_out = os.path.join(base_path, iid, "chr" + str(ch), prefix_out, "")
    if not os.path.exists(path_out):
            os.makedirs(path_out)
    ### Activate LOG FILE output if given
    if logfile == True:
        path_log = os.path.join(path_out, "hmm_run_log.txt")
        print(f"Set Output Log path: {path_log}")
        sys.stdout = open(path_log, 'w') 
    return path_out


def multi_run(fun, prms, processes = 4):
    """Implementation of running in Parallel.
    fun: Function
    prms: The Parameter Files
    processes: How many Processes to use"""
    print(f"Running {len(prms)} total jobs; {processes} in parallel.")
    
    with mp.Pool(processes = processes) as pool:
        results = pool.starmap(fun, prms)

def split_up_roh_df(base_path, path_out, iid, 
                    file_in="roh_info.csv", file_out="roh_gt.csv"):
    """Splits up the ROH-dataframe from base_path/file_in into file_out.
    Picks out Individual iid. Done to pass on "ground truth"
    base_path: Where to find roh_info.csv
    path_out: Where to save roh_gt to
    iid: Which Individual to extract from roh_info.csv."""
    #path = base_path + "roh_info.csv"
    path = os.path.join(base_path, file_in)
    dft = pd.read_csv(path, sep="\t")  # Load the Meta File

    save_df = dft[dft["iid"] == iid]
    save_path = os.path.join(path_out, file_out)
    save_df.to_csv(save_path, sep="\t", index=False)
    #print(f"Saved to {save_path}")
    return

def combine_individual_data(base_path, iid, delete=False, chs=range(1,23), 
                            prefix_out="", file_result="_roh_full.csv"):
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
    path_save = os.path.join(base_path, str(iid) + file_result)
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

######################################################
### For running bcftools & plink

def create_folders(input_base_folder, outfolder="plink_out/"):
    """Create Folders for ROH analysis with Plink/BCFTOOLs.
    Operates within HAPSBURG Mosaic Data Structure. Return
    h5 path, vcf path, and folder for intermediary output"""
    input_h5 = os.path.join(input_base_folder, "data.h5")
    input_vcf = os.path.join(input_base_folder, "data.vcf")
    
    if not os.path.exists(input_h5):
        raise RuntimeError(f"Create .vcf file: {input_h5}")
        
    plink_folder = os.path.join(input_base_folder, outfolder)
    if not os.path.exists(plink_folder):
        print(f"Creating Folder for: {plink_folder}")
        os.makedirs(plink_folder)
    
    return input_h5, input_vcf, plink_folder

def split_up_inferred_roh(df_t, iid, save_path):
    """Extract only ROH from Individual iid and saves it to save_path"""
    df_iid = df_t[df_t["iid"]==iid]
    df_iid.to_csv(save_path, index=False)
    
def postprocess_iid(df_plink, input_base_folder, iids, ch=3, prefix_out=""):
    """Split up results into roh.csv and roh_gt.csv for each IID.
    df_plink: Data Frame with Plink results, formated correctly"""

    for iid in iids:
        output_base_folder = os.path.join(input_base_folder, "output/")
        path_out = prepare_path(output_base_folder, iid, ch, prefix_out=prefix_out, logfile=False)

        path_inferred = os.path.join(path_out, "roh.csv")
        split_up_inferred_roh(df_plink, iid, save_path=path_inferred)   # Split up Inferred ROH
        split_up_roh_df(input_base_folder, path_out, iid)  # Split up Ground Truth ROH











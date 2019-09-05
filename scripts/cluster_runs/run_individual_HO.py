"""
Run Human Origin inference on the cluster.
Called with array jobs from sbatch
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import numpy as np
import os as os
import sys as sys
import pandas as pd
import socket

# Load the Meta File
meta_path = "../../Data/Marcus2019_1240k/meta_rev_unique_ids.csv"
meta_df = pd.read_csv(meta_path)
mod_df = meta_df[1098:]

sys.path.append("../../Python3/")  # Since now we are in the Root Directory
from hmm_inference import HMM_Analyze # Do not move. Should be after sys.path..

def prepare_path(path_output, iid, ch, prefix_out, logfile=True, output=False):
    """Prepare the path and pipe printing for one Individual
    logfile: Whether to pipe output to log-file"""
    # if not os.path.exists(path_output):
    #        raise RuntimeError(f"Path {path_output} not Found. Check!")
    path_log = os.path.join(path_output, str(
        iid), "chr" + str(ch), prefix_out, "")
    #path_log =  path_output + str(iid) + "/chr" + str(ch) + "/" + prefix_out

    if not os.path.exists(path_log):
        if output == True:
            print(f"Creating {path_log}...")
        os.makedirs(path_log)

    if logfile == True:
        path_log = path_log + "hmm_run_log.txt"
        if output == True:
            print(f"Set Output Log path: {path_log}")
        sys.stdout = open(path_log, 'w')


def analyze_chromosome_gt(iid, ch=3, n_ref=503, save=True, save_fp=False, exclude_pops=["TSI", ],
                          base_out_folder="./Empirical/HO/", prefix_out="gt/",
                          roh_in=100, roh_out=100, roh_jump=385, e_rate=0.01, e_rate_ref=0.001,
                          max_gap=0.0, logfile=True):
    """Run the analysis for one individual and chromosome on readcount data
    Wrapper for HMM Class. Takes 13 Parameters"""

    # The folder on what to run the Data on (Permanently set here to fixed loaction)
    h5_path_targets = "./Data/Marcus2019_1240k/mod_reich_sardinia_ancients_rev_mrg_dedup_3trm_anno.h5"
    # Path with the unique IDs per Modern Group
    meta_path_targets = "./Data/Marcus2019_1240k/meta_rev_unique_ids.csv"

    # Create Folder if needed, and pipe output if wanted
    prepare_path(base_out_folder, iid, ch, prefix_out, logfile=logfile)

    hmm = HMM_Analyze(cython=2, p_model="SardHDF5", e_model="diploid_gt", post_model="Standard",
                      manual_load=True, save=save, save_fp=save_fp)

    # Load and prepare the pre-processing Model
    hmm.load_preprocessing_model()              # Load the preprocessing Model
    hmm.p_obj.set_params(readcounts=False, destroy_phase=False,
                         prefix_out_data=prefix_out, excluded=exclude_pops, base_out_folder=base_out_folder,
                         h5_path_targets=h5_path_targets, meta_path_targets=meta_path_targets)

    # DELETE when run for with European Reference!!
    hmm.p_obj.set_params(h5_path1000g="./Data/1000Genomes/HDF5/1240kHDF5/all1240/chr",
                         meta_path_ref="./Data/1000Genomes/Individuals/meta_df_all.csv")

    hmm.load_data(iid=iid, ch=ch, n_ref=n_ref)  # Load the actual Data
    hmm.load_secondary_objects()

    # Set the Parameters
    hmm.e_obj.set_params(e_rate=e_rate, e_rate_ref=e_rate_ref)
    hmm.t_obj.set_params(roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump)
    hmm.post_obj.set_params(max_gap=max_gap)

    # hmm.calc_viterbi_path(save=save)           # Calculate the Viterbi Path.
    hmm.calc_posterior(save=save)              # Calculate the Posterior.
    hmm.post_processing(save=save)             # Do the Post-Processing.


#########################################################
def combine_individual_data(base_path, iid, delete=False, chs=range(1, 23), prefix_out=""):
    """Function to merge data from one Individual Analysis (all Chromosome)
    chs: Which Chromosomes to combine"
    delete: Whether to delete individual folder and contents after combining."""

    full_df_vec = []  # The full dataframe of inferred ROH blocks

    # Walk through Chromosomes and combine the Dataframes
    for ch in chs:
        path_roh = os.path.join(base_path, str(
            iid), "chr" + str(ch), prefix_out, "roh.csv")
        df_temp = pd.read_csv(path_roh, sep=",")
        full_df_vec.append(df_temp)

    full_df = pd.concat(full_df_vec)

    # Save to Path:
    path_save = os.path.join(base_path, str(iid) + "_roh_full.csv")
    full_df.to_csv(path_save, index=False)

    # Delete files in folder if need
    if delete == True:
        for ch in chs:
            path_folder = os.path.join(base_path, str(
                iid), "chr" + str(ch), prefix_out, "")

            for root, _, files in os.walk(path_folder):
                for file in files:
                    os.remove(os.path.join(root, file))
            os.rmdir(path_folder)  # Remove the Chromosome Folders
        # Remove the Individual Folder
        os.rmdir(os.path.join(base_path, str(iid), ""))

    return full_df

#########################################################


def analyze_individual_ho(iid, chs=range(1, 23), n_ref=2504, save=True, save_fp=False, exclude_pops=[],
                          base_out_folder="./Empirical/HO/", prefix_out="",
                          roh_in=100, roh_out=100, roh_jump=300, e_rate=0.001,
                          e_rate_ref=0.001, max_gap=0, logfile=True, output=True, delete=True):
    """Analyze a full single individual in a parallelized fasion. Run all Chromosome analyses in parallel
    Wrapper for analyze_chromosome_gt.
    logfile: Whether to use a logfile
    output: Whether to print general Output"""

    if output == True:
        print(f"Doing Individual {iid}...")

    # Prepare the Parameters for that Indivdiual
    prms = [[iid, ch, n_ref, save, save_fp, exclude_pops, base_out_folder, prefix_out,
             roh_in, roh_out, roh_jump, e_rate, e_rate_ref, max_gap, logfile] for ch in chs]

    # Run the analysis for all Parameters
    for prm in prms:
        analyze_chromosome_gt(*prms)

    # Merge results for that Individual
    combine_individual_data(base_out_folder, iid=iid, delete=delete, chs=chs)

    return

#########################################################
#########################################################

if __name__ == "__main__":
    analyze_individual_ho(iid="Koryak_5", chs=range(1,3), delete=False, logfile=False)

"""
Functions for analyzing Mosaic HDF5 Data
Contains Functions to run multiple Chromosomes and Individuals
"""

import pandas as pd
import numpy as np
import sys
sys.path.append("./Python3/")  # Since now we are in the Root Directory
from hmm_inference import HMM_Analyze   # Do not move. Should be after sys.path..

def analyze_individual(iid, ch=3, n_ref=503, save=True, save_fp=False,
                       path_mosaic="./Simulated/1000G_Mosaic/TSI/ch3_5cm/",
                       exclude_pops=["TSI", ], prefix_out=""):
    """Run the analysis for one individual and chromosome.
    Wrapper for HMM Class"""
    hmm = HMM_Analyze(cython=2, p_model="MosaicHDF5",
                      manual_load=True, save=save, save_fp=save_fp)

    # Load and prepare the pre-processing Model
    hmm.load_preprocessing_model()              # Load the preprocessing Model
    hmm.p_obj.set_folder(path_mosaic)         # Set the Folder
    hmm.p_obj.set_prefix_out_data(prefix_out)
    hmm.p_obj.set_exclude_pops(pops=exclude_pops)

    hmm.load_data(iid=iid, ch=ch, n_ref=n_ref)  # Load the actual Data
    hmm.load_emission_model()
    hmm.load_transition_model()

    hmm.set_diploid_observations()             # To diploidize Individuals
    hmm.t_obj.set_params(roh_in=1, roh_out=10, roh_jump=100)
    hmm.calc_viterbi_path(save=save)           # Calculate the Viterbi Path.
    hmm.calc_posterior(save=save)              # Calculate the Posterior.
    hmm.post_processing(save=save)             # Do the Post-Processing.
    print(f"Analysis of {iid} and Chr. {ch} successfully concluded!")


def run_full_individual(iid):
    """Run a Full Individual (all Chromosom).  NEEDS REWRITING"""
    ch_list = range(1, 23)
    for ch in ch_list:
        print(f"Doing Chromosome: {ch}")
        analyze_individual(iid=iid, ch=ch, save=True)


def split_up_df_roh_gt(base_path="../Simulated/1000G_Mosaic/TSI/ch3_10cm/", n=20, ch=3, prefix_out = ""):
    """Split up the Dataframe for ROH into the HMM Output Folders
    base_path: Outputfolder BaseException(
    ofile_name: The Name of the non-split dataframe at base
    n: How many individuals
    prefix_out: String that gets inserted before file in save path
    """

    path = base_path + "roh_info.csv"
    dft = pd.read_csv(path, sep="\t")  # Load the Meta File

    iids = ["iid" + str(i) for i in range(n)]

    for iid in iids:
        save_df = dft[dft["iid"] == iid]
        save_path = base_path + "output/" + \
            iid + "/chr" + str(ch) + "/" + prefix_out + "roh_gt.csv"
        save_df.to_csv(save_path, sep="\t", index=False)


def multi_run_individuals(n=20, ch=3, base_path="./Simulated/1000G_Mosaic/TSI/",
                          lengths=[1, 3, 5, 10], exclude_pops=["TSI", ], prefix_out=""):
    """Code for analyzing multiple Simulated Indiviuals.
    n: How many per run
    ch: Which chromosome
    base_path: Which Path (Including Population)
    lengths: Which Lengths (0 is no block).
    exclude_pops: Which populations to exclude as reference"""

    iids = ["iid" + str(i) for i in range(n)]
    folders = [base_path + "ch"
               + str(ch) + "_" + str(int(l)) + "cm/" for l in lengths]

    for f in folders:
        for iid in iids:
            analyze_individual(iid=iid, path_mosaic=f,
                               exclude_pops=["TSI", ], prefix_out=prefix_out)

        # Could do before for full parallelization
        #split_up_df_roh_gt(base_path=f, n=n, ch=ch, prefix_out=prefix_out)


if __name__ == "__main__":
    #path_mosaic = "./Simulated/1000G_Mosaic/TSI/ch3_10cm/"
    #analyze_individual(iid="iid0", path_mosaic=path_mosaic)

    # multi_run_individuals(n=20, ch=3, base_path="./Simulated/1000G_Mosaic/TSI/",
    #                      lengths=[1, 3, 5, 10])

    multi_run_individuals(n=20, ch=3, base_path="./Simulated/1000G_Mosaic/TSI/",
                          lengths=[1, 3, 5, 10, ], prefix_out= "baseline/")

    #multi_run_individuals(n=100, ch=3, base_path="./Simulated/1000G_Mosaic/TSI/", lengths=[0,])

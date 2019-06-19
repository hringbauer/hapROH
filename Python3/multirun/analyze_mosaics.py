"""
Functions for analyzing Mosaic HDF5 Data
Contains Functions to run multiple Chromosomes and Individuals
"""
import pandas as pd
import numpy as np
import sys
sys.path.append("./Python3/")  # Since now we are in the Root Directory
from hmm_inference import HMM_Analyze  # Do not move


def analyze_individual(iid, ch=3, n_ref=503, save=True, save_fp=False, path_mosaic="./Simulated/1000G_Mosaic/TSI/ch3_5cm/"):
    """Run the analysis for one individual and chromosome.
    Wrapper for HMM Class"""
    hmm = HMM_Analyze(cython=2, p_model="MosaicHDF5",
                      manual_load=True, save=save, save_fp=save_fp)

    hmm.load_preprocessing_model()              # Load the preprocessing Model
    hmm.p_obj.set_folder(path_mosaic)         # Set the Folder
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

def multi_run_individuals():
    """Code for running multiple Indiviuals (Folders)"""
    n=20 # Nr of analyzed Individuals
    ch=3
    iids = ["iid" + str(i) for i in range(n)]
    lengths = [1, 3, 5, 10]
    base_path = "./Simulated/1000G_Mosaic/TSI/"  # Start of SavePaths

    folders = [base_path + "ch" + str(ch) + "_" + str(int(l)) + "cm/" for l in lengths]

    for f in folders:
        for iid in iids:
            analyze_individual(iid=iid, path_mosaic=f)

if __name__ == "__main__":
    path_mosaic = "./Simulated/1000G_Mosaic/TSI/ch3_10cm/"
    #analyze_individual(iid="iid0", path_mosaic=path_mosaic)
    multi_run_individuals()

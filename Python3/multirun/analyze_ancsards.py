"""
Functions for analyzing the ancient Sardinian Data
Contains Functions to run multiple Chromosomes and Individuals
"""

import sys
sys.path.append("./Python3/") # Since now we are in the Root Directory

from hmm_inference import HMM_Analyze
import numpy as np
import pandas as pd

############## Useful Parameters

meta_path = "./../ancient-sardinia/output/meta/meta_final.csv"
anc_ind = 1029  # Up to which individual there are ancestrals
anc_sardind = 57  # Nr of ancient Sarinian individuals

############# The Functions

def analyze_individual(iid, ch, n_ref=503, save=True, save_fp=False):
    """Run the analysis for one individual and chromosome.
    Wrapper for HMM Class"""
    hmm = HMM_Analyze(cython=2, p_model="SardHDF5",
                      manual_load=True, save=save, save_fp=save_fp)

    hmm.load_objects(iid=iid, ch=ch, n_ref=n_ref)
    hmm.set_diploid_observations()
    hmm.t_obj.set_params(roh_in=1, roh_out=10, roh_jump=100)
    hmm.calc_viterbi_path(save=save)           # Calculate the Viterbi Path.
    hmm.calc_posterior(save=save)              # Calculate the Posterior.
    hmm.post_processing(save=save)             # Do the Post-Processing.
    print(f"Analysis of {iid} and Chr. {ch} successfully concluded!")


def run_full_individual(iid):
    """Run a Full Individual (all Chromosom)"""
    ch_list = range(1, 23)
    for ch in ch_list:
        print(f"Doing Chromosome: {ch}")
        analyze_individual(iid=iid, ch=ch, save=True)


def run_all_sardinians(min_cov=0.5):
    """Run all Sardinian Ancients"""
    meta_df = pd.read_csv(meta_path)  # Load the Meta File
    as_df = meta_df[:anc_sardind]
    iid_list = as_df[(as_df["mean_cov"] > min_cov)]["iid"].values
    print(f"Found n={len(iid_list)} Inds with Cov > {min_cov}")

    for iid in iid_list:
        print(f"Doing iid: {iid}")
        run_full_individual(iid)

    print("\nFinished all. GZ")


if __name__ == "__main__":
    analyze_individual(iid="MA89", ch=6, save=True, save_fp=False)
    #run_all_sardinians()

from hmm_inference import analyze_individual
import numpy as np
import pandas as pd

meta_path = "./../ancient-sardinia/output/meta/meta_final.csv"
anc_ind = 1029  # Up to which individual there are ancestrals
anc_sardind = 57  # Nr of ancient Sarinian individuals

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
    run_all_sardinians()

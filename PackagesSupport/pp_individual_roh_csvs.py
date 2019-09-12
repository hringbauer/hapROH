"""
Helper Functions to post-process large Number of Individuals on Cluster into Summary ROH data tables
(saved as pandas)
Structure: Produce list of paths with Individuals ROH files for a group.
These are then combined into one summary data frame (and some post-processing such as gap merging can be done)
This is then saved as a single "summary" dataframe (csv) 
@ Author: Harald Ringbauer, 2019, All rights reserved
"""
import numpy as np
import pandas as pd
import os


def give_iid_paths(iids, base_folder="./Empirical/HO/", suffix = "_roh_full.csv"):
    """Return list of the paths to each ROH.csv.
    Combine basefolder and iids."""
    paths_csv = [os.path.join(base_folder, str(iid) + suffix) for iid in iids]
    return paths_csv

def create_combined_ROH_df(paths, iids, pops, min_cm=4, snp_cm=100, 
                     savepath = "./Empirical/HO/roh_summary_HO.csv", gap=0.5, output=True):
    """Create the Ancient Sardinian Summary Dataframe
    paths: List of .csv Paths to load the Data from
    snp_cm: Minimum SNP Density per cM
    savepath: If given, where to save the summary .csv to
    gap: Gaps to merge [in cM]"""
    n = len(paths)

    # Default Values: Zero ROH
    max_roh = np.zeros(n, dtype="float")
    sum_roh = np.zeros(n, dtype="float")
    n_roh = np.zeros(n, dtype="int")

    for i, path in enumerate(paths):
        df_roh = pd.read_csv(path)
        if gap>0:
            df_roh = merge_called_blocks(df_roh, max_gap=gap/100)
        df_roh = post_process_roh_df(df_roh, output=output, snp_cm=snp_cm, min_cm=min_cm)  ### Do the Postprocessing
        max_roh[i], sum_roh[i], n_roh[i] = individual_roh_statistic(df_roh, output=output)

    ### Create the Dataframe:
    d = {"iid": iids, "pop" : pops, "max_roh": max_roh*100, "sum_roh" : sum_roh*100, "n_roh" : n_roh}
    df1 = pd.DataFrame(d).sort_values(by="sum_roh", ascending=False)  # Sort output
    
    if len(savepath) > 0:
        df1.to_csv(savepath, sep=",", index=False)
        print(f"Saved to: {savepath}")
        
    return df1

def merge_called_blocks(df, max_gap=0, output=False):
        """Merge Blocks in Dataframe df and return merged Dataframe"""
        if len(df) == 0:
            return df  # In case of empty dataframe don't do anything

        df_n = df.drop(df.index)  # Create New Data frame with all raws removed
        row_c = df.iloc[0, :].copy()

        # Iterate over all rows, update blocks if gaps small enough
        for index, row in df.iterrows():
            if (row["StartM"] - row_c["EndM"] < max_gap) and (row["ch"] == row_c["ch"]):
                row_c["End"] = row["End"]
                row_c["EndM"] = row["EndM"]
                row_c["length"] = row_c["End"] - row_c["Start"]
                row_c["lengthM"] = row_c["EndM"] - row_c["StartM"]

            else:  # Save and go to next row
                df_n.loc[len(df_n)] = row_c  # Append a row to new df
                row_c = row.copy()

        df_n.loc[len(df_n)] = row_c   # Append the last row

        if output == True:
            print(f"Merged n={len(df) - len(df_n)} gaps < {max_gap} M")
        return df_n








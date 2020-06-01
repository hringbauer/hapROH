"""
Helper Functions to post-process large Number of Individuals on Cluster into Summary ROH data tables
(saved as pandas)
Structure: Produce list of paths with Individuals ROH files for a group.
These are then combined into one summary data frame (and some post-processing such as gap merging can be done)
This is then saved as a single "summary" dataframe (csv) 
@ Author: Harald Ringbauer, 2019
"""
import numpy as np
import pandas as pd
import os

def give_iid_paths(iids, base_folder="./Empirical/HO/", suffix = "_roh_full.csv"):
    """Return list of the paths to each ROH.csv.
    Combine basefolder and iids."""
    paths_csv = [os.path.join(base_folder, str(iid) + suffix) for iid in iids]
    return paths_csv

def create_combined_ROH_df(paths, iids, pops, min_cm=[4,8,12], snp_cm=100, 
                           gap=0.5, min_len1=2, min_len2=4, output=True, sort=True):
    """Create the Ancient Sardinian Summary Dataframe
    paths: List of .csv Paths to load the Data from
    snp_cm: Minimum SNP Density per cM
    min_cm: Minimal centiMorgan for Postprocessing Postanalysis
    savepath: If given, where to save the summary .csv to
    gap: Maximum Gaps to merge [in cM]
    min_len1: Minimum Length of both Blocks for merge [in cM]
    min_len2: Minimum Length of longer Block to merge [in cM]"""    
    
    ### Dry run to check if paths exist
    idcs, paths_exist, paths_not_exist = [], [], []
    for i, path in enumerate(paths):
        if os.path.exists(path):
            paths_exist.append(path)
            idcs.append(i)
        else:
            paths_not_exist.append(path)
            
    if len(paths_not_exist)>0:
        print(f"Warning, could not find {len(paths_not_exist)} Paths:")
        print(paths_not_exist)
    
    df_rohs = [pd.read_csv(path) for path in paths_exist] # Create list of Dataframes
    # A bit wast of memory, but allows easy refactoring
    iids, pops = iids[idcs], pops[idcs]  # Cut out only the existing Indivdiduals
    
    df1 =combine_ROH_df(df_rohs, iids=iids, pops=pops, min_cm=min_cm, snp_cm=snp_cm,
                        gap=gap, min_len1=min_len1, min_len2=min_len2,
                        output=output, sort=sort)
    return df1


def combine_ROH_df(df_rohs, iids=[], pops=[], min_cm=[4,8,12], snp_cm=100, 
                   gap=0.5, min_len1=2, min_len2=4, output=True, sort=True):
    """Takes list of ROH Dataframes, and creates a single
    summary dataframe. Return single dataframe.
    Being wrapped around by create_combined_ROH_df
    which does the path and loading.
    df_rohs: List of individual dataframes
    iids: IIds of the Individuals (filled in columns)
    pops: Populations of the Individuals (filled in column)
    gap, min_len1/2 are in cM (!!)"""
                   
    # Create Place Holder Arrays
    n = len(df_rohs)
    l = len(min_cm)
    max_roh = np.zeros(n, dtype="float")
    sum_roh = np.zeros((n,l), dtype="float")
    n_roh = np.zeros((n,l), dtype="int")
                        
    for i, df_roh in enumerate(df_rohs):
        if gap>0:
            df_roh = merge_called_blocks(df_roh, max_gap=gap/100, 
                                         min_len1=min_len1/100, min_len2=min_len2/100)

        for j, m_cm in enumerate(min_cm):
            df_roh = post_process_roh_df(df_roh, output=output, snp_cm=snp_cm, min_cm=m_cm)  ### Do the Postprocessing
            max_roh_t, sum_roh[i,j], n_roh[i,j] = individual_roh_statistic(df_roh, output=output)
            if j==0:  # Only calcuate maximum for shortest cM category (containing the bigger ones)
                max_roh[i] = max_roh_t

    ### Create the Dataframe with the basic entries:
    d = {"iid": iids, "pop" : pops, "max_roh": max_roh*100}
    df1 = pd.DataFrame(d) 
    
    ### Add Values for varying cM cutoffs:
    for j, m_cm in enumerate(min_cm):
        df1[f"sum_roh>{m_cm}"]= sum_roh[:,j]*100
        df1[f"n_roh>{m_cm}"]=n_roh[:,j]
    if sort:
        df1 = df1.sort_values(by=f"sum_roh>{min_cm[0]}", ascending=False)  # Sort output by minimal Cutoff                    
    return df1                 
                      
def merge_called_blocks(df, max_gap=0, min_len1=0.02, 
                        min_len2=0.04, output=False):
        """Merge Blocks in Dataframe df and return merged Dataframe.
        Gap is given in Morgan"""
        if len(df) == 0:
            return df  # In case of empty dataframe don't do anything

        df_n = df.drop(df.index)  # Create New Data frame with all raws removed
        row_c = df.iloc[0, :].copy()
        #row_c["lengthM"] = row_c["EndM"] - row_c["StartM"] # Should be there

        # Iterate over all further rows, update blocks if gaps small enough
        for index, row in df.iloc[1:,:].iterrows():
            ### Calculate Conditions
            short_b = np.min([row_c["lengthM"], row["lengthM"]])>=min_len1
            long_b = np.max([row_c["lengthM"], row["lengthM"]])>=min_len2
            short_g = (row["StartM"] - row_c["EndM"] < max_gap)
            same_c = (row["ch"] == row_c["ch"])
            
            if short_b and long_b and short_g and same_c:
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
    
def post_process_roh_df(df, min_cm=4, snp_cm=60, output=True):
    """Post Process ROH Dataframe. Filter to rows that are okay.
    min_cm: Minimum Length in CentiMorgan
    snp_cm: How many SNPs per CentiMorgan"""
    densities = df["length"] / (df["lengthM"] * 100)
    densities_ok = (densities > snp_cm)
    
    df["SNP_Dens"] = densities
    
    # Filter for SNP Density:
    df = df[densities_ok]
    
    # Filter for Length:
    length_okay = (df["lengthM"] * 100) > min_cm
    
    if output==True:
        print(f"Min SNPs per cM> {snp_cm}: {np.sum(densities_ok)}/{len(densities_ok)}")
        print(f"> {min_cm} cM: {np.sum(length_okay)}/{len(length_okay)}")
    
    df = df[length_okay]
    return df

def individual_roh_statistic(df, output=True):
    """Gives out Summary statistic of ROH df"""
    if len(df)==0:   # If no ROH Block found
        max_roh, sum_roh, n_roh = 0, 0, 0
    else:
        max_roh = np.max(df["lengthM"])
        sum_roh = np.sum(df["lengthM"])
        n_roh = len(df["lengthM"])
     
    if output==True:
        print(f"Max. ROH: {max_roh * 100:.3f}")
        print(f"Sum. ROH: {sum_roh * 100:.3f}")
        print(f"Nr. ROH: {n_roh}")   
    return max_roh, sum_roh, n_roh

#########################################################
#### Main Function to post-process Individual ROH outputs

def pp_individual_roh(iids, meta_path="./Data/ReichLabEigenstrat/Raw/meta.csv", 
                      base_folder="./Empirical/Eigenstrat/Reichall/", 
                      suffix='_roh_full.csv', save_path="", min_cm=[4,8,12], snp_cm=50, 
                      gap=0.5, min_len1=2.0, min_len2=4.0,
                      output=True, meta_info=True):
    """Post-process Individual ROH .csv files. Combines them into one summary ROH.csv, saved in save_path.
    Use Individuals iids, create paths and run the combining.
    iids: List of target Individuals
    base_folder: Folder where to find individual results .csvs
    min_cm: Minimum post-processed Length of ROH blocks. Array (to have multiple possible values)
    snp_cm: Minimum Number of SNPs per cM
    gap: Maximum length of gaps to merge
    output: Whether to plot output per Individual.
    meta_info: Whether to merge in Meta-Info from the original Meta File
    save_path: If given, save resulting dataframe there
    min_len1: Minimum Length of shorter block to merge [cM]
    min_len2: Maximum Length of longer block to merge [cM]"""
    
    ### Look up Individuals in meta_df and extract relevant sub-table
    df_full = pd.read_csv(meta_path)
    df_meta = df_full[df_full["iid"].isin(iids)]  # Extract only relevant Indivdiuals
    
    print(f"Loaded {len(df_meta)} / {len(df_full)} Individuals from Meta")
    
    paths = give_iid_paths(df_meta["iid"], base_folder=base_folder, suffix=suffix)
    df1 = create_combined_ROH_df(paths, df_meta["iid"].values, df_meta['clst'].values, 
                                 min_cm=min_cm, snp_cm=snp_cm, gap=gap, 
                                 min_len1=min_len1, min_len2=min_len2, output=output)
    
    ### Merge results with Meta-Dataframe
    if meta_info:
        df1 = pd.merge(df1, df_meta, on="iid")
        
    if len(save_path) > 0:
        df1.to_csv(save_path, sep="\t", index=False)
        print(f"Saved to: {save_path}")
    
    return df1

##########################################################
### Postprocess Dataframes based on Geography and Keywords

####################################################
##### Helper Functions to split up Dataframes

def extract_df_geo(df, lat0, lat1, lon0, lon1):
    """Extract Dataframe df from Sub Data frame based on coordinates
    lat0,lat1: Min and Max Lat. Equ. for lon0,lon1"""
    lat_okay = (df["lat"]>lat0) & (df["lat"]<lat1)
    lon_okay = (df["lon"]>lon0) & (df["lon"]<lon1)
    df_s = df[lat_okay & lon_okay]
    return df_s

def extract_df_age(df, age0, age1=1e6):
    """Extract Dataframe based on age.
    df: Input Dataframe; age0 and age1 min and max age"""
    age_okay = (df["age"]>=age0) & (df["age"]<=age1)
    df = df[age_okay]
    return df

def give_df_clsts(df, search=[], col="pop"):
    """Return sub dataframe within df
    where in col one of search strings (list of string)"""
    if len(search)>0:
        idx = df[col].str.contains('|'.join(search)) # Find
        df = df[idx]
    else:
        df=df[0:0]
    return df

def extract_sub_df_geo_kw(df, lat0, lat1, lon0, lon1, keywords=[], output=True):
    """Extract Dataframe df from Sub Data frame based on coordinates
    AND from keywords
    lat0,lat1: Min and Max Lat. Equ. for lon0,lon1"""
    df1 = extract_df_geo(df, lat0, lat1, lon0, lon1).copy()  # Make a new Copy. America
    df2 = give_df_clsts(df, search=keywords) 
    df_m = pd.concat([df1, df2]).drop_duplicates().reset_index(drop=True)
    if output:
        print(f"Found {len(df_m)} Individuals; {len(df1)} from Geography")
    return df_m


#################################################################
##### More utility functions

def calc_average_roh(df, cms=[4,8,12,20], 
                     col_pop="pop", new_pop="Average"):
    """Calcualte the average ROH in df.
    in columns cmss [list]
    Return new dataframe"""
    cs = [4,8,12,20]
    labels = [f"sum_roh>{c}" for c in cs]

    df_new = pd.DataFrame({"iid": ["Average"], col_pop:[new_pop]})
    for l in labels:
        df_new[l] = np.mean(df[l])
    return df_new







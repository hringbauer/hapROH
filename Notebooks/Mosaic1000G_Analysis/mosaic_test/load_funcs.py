"""
Classes and functions to Load the simulated and processed ROH Data
Python Functions to keep functions and analysis seperate
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import pandas as pd
import numpy as np

class Summary_ROH_Calls(object):
    """Class to Provide Summary statistics for multiple replicate Runs"""
    mosaic_folder = "../../Simulated/1000G_Mosaic/TSI/"
    output_prefix = ""
    ch = 3 ### Which Chromosome to 
    nr_iid = 20 # How many Individuals were simulated
    blen_cm = 1 # The Length of the Block to analyze next
    error = 0.0   # For the Subclasses
    missing = 0 # For the Subclasses
    output = 1 # The Level of output. 0 None, 1 All,

    def __init__(self, mosaic_folder, ch=3, nr_iid = 20, blen_cm = 1, output_prefix = "", 
                 error=0, missing=0, output=0):
        """Initialize the whole Class"""
        self.mosaic_folder = mosaic_folder
        self.ch = ch
        self.blen_cm = blen_cm
        self.nr_iid = nr_iid
        self.output_prefix = output_prefix
        self.output=output
        self.error=error
        self.missing=missing

    def provide_iid_folders(self):
        """Return a list of folders into which replicate simulations are saved into
        folder: Full folder of the analysis
        ch: Which Chromosome was analyzed
        nr_iid: How many replicate Individual were produced
        bl_cm: The length of the analyzed block class in cM
        output_prefix: Prefix before the File
        mode: Which Folder Mode to load. 0: Standard 1: Error Folders 2: Missing Folders"""
        folder, ch, nr_iid, blen_cm, output_prefix = self.mosaic_folder, self.ch,  self.nr_iid, self.blen_cm, self.output_prefix
        
        path_1 = "ch" + str(ch) + "_" + str(int(blen_cm)) + "cm/output/"  # The first part of the Path
        path_2 = "/chr" + str(ch) + "/" + output_prefix # The last part of the Path (within the folder)
        iid_list = ["iid" + str(i) for i in range(nr_iid)]  # The middle parts of the path

        full_paths = [(folder + path_1 + str(i) + path_2)  for i in iid_list]
        return full_paths
    
    def collect_power_df(self):
        """Create and return the Dataframe with the Power.
        Check every simulated block for overlap."""
                
        folders = self.provide_iid_folders()   # Load all the folders
        df_calls = []

        n_call, n_sim = 0, 0  # Number of total called blocks and simulated blocks
        for f in folders:   
            df_o = load_observed(f)   # Loading should throw error if not existend
            df_s = load_simulated(f)

            n_call += len(df_o)
            n_sim += len(df_s)

            df_call = check_obs_vrs_sim(df_o, df_s)
            df_calls.append(df_call)

        df_calls = pd.concat(df_calls)   # Concatenate all the Results

        if self.output==1:
            print(f"Total Nr Simulated: {n_sim}")
            print(f"Total Nr Called: {n_call}")
            
        return df_calls

    def collect_fp_df(self):
        """Collect and return the Dataframe with the false positive Calls"""
        folders = self.provide_iid_folders()   # Load all the folders
        observed_dfs = [load_observed(f) for f in folders]
        df_observed = pd.concat(observed_dfs)
        return df_observed
    
    #########################################
    ### Give back arrays of information, needed for plots and analysis of multiple folders
    
    def give_power_dfs(self, bl_lens, ovlp_frac=0.8):
        """Load the Power dfs for vector of Block Lengths bl_lens
        bl_lens: Array with length of Blocks
        Return list of dataframes of all blocks, list of dataframe of all called blocks, and list of power"""
        mosaic_folder, output_prefix = self.mosaic_folder, self.output_prefix
        
        df_call_vec = []
        
        for l in bl_lens:
            self.blen_cm = l
            df_call_vec.append(self.collect_power_df())
            
        # Post -Process
        powers = [calc_power(df, ovlp_frac=ovlp_frac) for df in df_call_vec]   # Calculate the Power per Block Length
        df_called = [return_calls_only(df) for df in df_call_vec]   # Only Keep the called blocks for plotting
        return df_call_vec, df_called, powers
    
#############################################
    
class Summary_ROH_Calls_Error(Summary_ROH_Calls): 
    """Same as Summary Calls but with updated folder function"""
    error = 0
    
    def provide_iid_folders(self):
        """Return a list of folders into which replicate simulations are saved intos"""
        folder, ch, nr_iid, blen_cm, output_prefix = self.mosaic_folder, self.ch,  self.nr_iid, self.blen_cm, self.output_prefix
        error = self.error

        e_print = str(round(error, 4)).split(".")[1] # Extract four digits after decimal 
        path_1 = "ch" + str(ch) + "_" + str(int(blen_cm)) + "cm/error/" + e_print + "/output/"  # The first part of the Path
        path_2 = "/chr" + str(ch) + "/" + output_prefix # The last part of the Path (within the folder)
        iid_list = ["iid" + str(i) for i in range(nr_iid)]  # The middle parts of the path

        full_paths = [(folder + path_1 + str(i) + path_2)  for i in iid_list]
        return full_paths
    
#############################################
    
class Summary_ROH_Calls_Missing(Summary_ROH_Calls):
    """Same as Summary Calls but with updated folder function"""
    
    def provide_iid_folders(self):
        """Return a list of folders into which replicate simulations are saved intos"""
        folder, ch, nr_iid, blen_cm, output_prefix = self.mosaic_folder, self.ch,  self.nr_iid, self.blen_cm, self.output_prefix
        missing = self.missing
        m_print = str(round(missing, 4)).split(".")[1] # Extract four digits after decimal 
        path_1 = "ch" + str(ch) + "_" + str(int(blen_cm)) + "cm/missing/" + m_print + "/output/"  # The first part of the Path
        path_2 = "/chr" + str(ch) + "/" + output_prefix # The last part of the Path (within the folder)
        iid_list = ["iid" + str(i) for i in range(nr_iid)]  # The middle parts of the path
        full_paths = [(folder + path_1 + str(i) + path_2)  for i in iid_list]
        return full_paths

#############################################
#############################################
### Various Helper Functions for Analysis

def load_observed(path, file="roh.csv"):
    """ Load simulated Dataframe from path"""
    df = pd.read_csv(path + file, sep=",")  # Load the Meta File
    return df

def load_simulated(path, file="roh_gt.csv"):
    """ Load Ground Truth Dataframe from path"""
    df = pd.read_csv(path + file, sep="\t")  # Load the Meta File
    return df

def find_overlap(l, min_l, max_l):
    """ Find overlap of Interval l with intervals starting at min_l and ending at max_l
    Return max overlap as well as Block length (or 0, 0 if no overlap)
    l: Interval [Length 2 list]
    min_l: Minimum Interval Lengths
    max_l: Maximum Interval Lengths"""
    assert(len(l)==2)
    assert(len(min_l) == len(max_l))

    if len(min_l)==0: # If no interval given return 0 Overlap 0 total Length
        return 0, 0

    min_both = np.maximum(l[0], min_l)
    max_both = np.minimum(l[1], max_l)  
    overlap = max_both - min_both   # Calculate the Overlap

    i = np.argmax(overlap)
    max_overlap = overlap[i]

    if max_overlap < 0:
        return 0., 0.   # Return 0 Overlap and 0 total length

    else:
        orginal_length = max_l[i] - min_l[i]
        return max_overlap, orginal_length


def check_obs_vrs_sim(df_o, df_s):
    """Analyze observed . 
    Input: 2 Dataframes
    Output: 1 Dataframe
    Return: One dateframe with intersecting [Overlap, CalledLength, Originallength, Position]"""
    begin_obs = df_o["StartM"].values * 100  # Factor 100: To do everything in centiMorgan
    end_obs = df_o["EndM"].values * 100

    overlaps = np.zeros(len(df_s))
    lengths = np.zeros(len(df_s))

    for i, row in df_s.iterrows():
        roh = [row["ROH_Begin"] * 100, row["ROH_End"] * 100]
        ovlp, lgth = find_overlap(roh, begin_obs, end_obs)  # Find the Overlap

        overlaps[i] = ovlp
        lengths[i] = lgth

    org_lengths = (df_s["ROH_End"] - df_s["ROH_Begin"]) * 100 ### The original Lengths
    ovlp_frac = overlaps / org_lengths ### Fraction of Block called

    df_call = pd.DataFrame({"Overlap": overlaps, "CalledLength":lengths, "OriginalLength": org_lengths,
                            "OverlapFrac": ovlp_frac, "Position": df_s["ROH_Begin"].values * 100, "iid":df_s["iid"]})

    ###  Create df with Original Length, Found Length, Position, Overlap (Fraction)  
    return df_call

#############################################
### Post-Processing the Call Dataframe

def statistics_power_df(df, min_frac=0.8):
    """Report some statistics on the Power Dataframe.
    min_len: Minimum Overlap Length of Block to be counted as a call"""
    l = len(set(df["iid"]))
    no_calls = df["Overlap"]<0.01
    df_c = df[~no_calls]  # The Dataframe with the lenghts called
    mean_called = np.mean(df_c["CalledLength"])
    good_call_nr = np.sum(df["OverlapFrac"] > min_frac)

    print(f"{l} unique Individuals")
    print(f"{len(df) - np.sum(no_calls)} / {len(df)} Blocks called")
    print(f"{good_call_nr} Blocks called > {min_frac*100} %")
    print(f"{mean_called:.4f} average Block Length cM (called)")

def calc_power(df, ovlp_frac=0.8):
    """Calculate and return the Power from the Calling Dataframe"""
    no_calls = np.sum(df["OverlapFrac"] < ovlp_frac)
    power = 1  - (no_calls / len(df))
    return power

def give_SE(df):
    """Calculate and return Standard Deviation (and maybe quartiles) [in cM]"""
    s = np.std(df["CalledLength"].values)
    return s

def give_bias(df, true_lgth):
    """Calculate and return the bias of estimates.
    true_length: Length of simulated blocks [in cM]"""
    b = np.mean(df["CalledLength"].values) - true_lgth
    return b

def return_calls_only(df, ovlp=0.01):
    """Return Dataframe of only called Blocks"""
    no_calls = df["Overlap"]<ovlp
    df_c = df[~no_calls]  # The Dataframe with the lenghts called
    return(df_c)

def false_power_statistics(df):
    '''Report some statistics on the false Power Dataframe'''
    bl_lens = df["lengthM"].values

    print(f"Found {len(bl_lens)} FP blocks")
    if len (bl_lens)>0:
        print(f"Average Block length: {np.mean(bl_lens):.4f} cM")
        print(f"Maximum Block length: {np.max(bl_lens):.4f} cM")
        
def false_positive_nrs(df_call_fp, pw_lens):
    """For list pw_lens (in cM) return the number of blocks longer than each of that
    Return array of integers, same in size as pw_lens"""
    lengths = df_call_fp["lengthM"] * 100
    powers = [np.sum(lengths>l) for l in pw_lens]
    return np.array(powers) 

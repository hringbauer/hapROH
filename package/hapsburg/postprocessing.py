"""
Class for calling ROH from Posterior Data. Saves results as a .csv
created by Pandas.
Contains Sub-Classes, as well as factory Method.
Pls always change parameters with set_params method!
@ Author: Harald Ringbauer, 2019
"""

import numpy as np
import pandas as pd
import os


class PostProcessing(object):
    """Class that does PostProcessing of HAPSBURG output.
    Has Methods to save the output as """
    folder = ""          # The Folder to operate in
    roh_df = []          # Dataframe for Runs of Homozygosity

    cutoff_post = 0.8  # Cutoff Probability for ROH State
    roh_min_l = 0.01  # Cutoff [in Morgan]
    max_gap = 0.01  # The Maximum Gap Length to be Merged [in Morgan]
    snps_extend = 0  # SNPs to snip on either side of blocks (after merge)

    merge = True  # Whether to Merge ROH Blocks
    output = True
    save = True  # Whether to save output into Folder

    def __init__(self, folder="", load=False, output=True, save=True):
        """Initialize Class.
        Load: Whether to immediately Load the Posterior Data"""
        self.folder = folder
        self.output = output
        self.save = save

        if load == True:
            self.load_data()

    def set_params(self, **kwargs):
        """Set the Parameters.
        Takes keyworded arguments"""
        for key, value in kwargs.items():
            setattr(self, key, value)

    def load_data(self, folder=""):
        """Load and return genetic Map [l], 
        positions [l] and Posterior0 [l]"""
        if len(folder) == 0:
            folder = self.folder  # Use the Folder of the Class

        # Load Posterior
        post_path = folder + "posterior0.csv"
        posterior0 = np.loadtxt(post_path, dtype="float", delimiter=",")

        # Load Linkage Map
        map_path = folder + "map.csv"
        pos_path = folder + "pos.csv"

        if os.path.exists(map_path):
            r_map = np.loadtxt(
                map_path, dtype="float", delimiter=",")
        else:
            print("No Genetic Map found!!! Defaulting...")
            r_map = np.arange(len(self.posterior0))
        
        if os.path.exists(pos_path):
            pos = np.loadtxt(
                pos_path, dtype="int", delimiter=",")
        else:
            print("No physical Positions found!!! Defaulting...")
            pos = np.arange(len(self.posterior0))

        assert(len(r_map) == len(posterior0))  # Sanity Check
        assert(len(pos) == len(posterior0))  # Sanity Check
        print(f"Successfully loaded for PP. from {folder}")

        return r_map, pos, posterior0

    def merge_called_blocks(self, df, max_gap=0):
        """Merge Blocks in Dataframe df and return merged Dataframe"""
        if len(df) == 0:
            return df  # In case of empty dataframe don't do anything

        if max_gap == 0:
            max_gap = self.max_gap

        df_n = df.drop(df.index)  # Create New Data frame with all raws removed
        row_c = df.iloc[0, :].copy()

        # Iterate over all rows, update blocks if gaps small enough
        for index, row in df.iterrows():
            if row["StartM"] - row_c["EndM"] < max_gap:
                row_c["End"] = row["End"]
                row_c["EndM"] = row["EndM"]
                row_c["length"] = row_c["End"] - row_c["Start"]
                row_c["lengthM"] = row_c["EndM"] - row_c["StartM"]

            else:  # Save and go to next row
                df_n.loc[len(df_n)] = row_c  # Append a row to new df
                row_c = row.copy()

        df_n.loc[len(df_n)] = row_c   # Append the last row

        if self.output == True:
            print(f"Merged n={len(df) - len(df_n)} gaps < {max_gap} M")
        return df_n

    def snp_extend(self, df, r_map, snps_extend=0):
        """Extend Blocks in df by # n_snps.
        Return dataframe of same size but with modified blocks"""
        if snps_extend == 0:
            snps_extend = self.snps_extend   # Default with Object Value

        starts = df["Start"].values - snps_extend
        ends = df["End"].values + snps_extend
        starts = np.maximum(0, starts)  # To not shoot out of boundary
        ends = np.minimum(len(r_map), ends)  # Same
        ends = np.maximum(starts, ends)  # Make sure that ends are after starts
        l = ends - starts

        ends_map = r_map[ends - 1]  # -1 to stay within bounds
        starts_map = r_map[starts]
        l_map = ends_map - starts_map
        l_diff = l_map - df["lengthM"].values

        # Set the new Values
        df["Start"] = starts
        df["End"] = ends
        df["StartM"] = starts_map
        df["EndM"] = ends_map
        df["lengthM"] = l_map
        df["length"] = l

        if self.output == True:
            print(f"Extended n={len(df)} ROH by {snps_extend} SNPs on Avg. by {np.mean(l_diff):.6f} M")
        return df

    def modify_posterior0(self, posterior0):
        """Load and return the posterior."""
        #roh_post = 1 - np.exp(posterior0)  # Go to non-logspace probability
        roh_post = 1 - posterior0  # Go to non-logspace probability
        return roh_post
    
    def create_df(self, starts, ends, starts_map, ends_map, 
                  l, l_map, iid, ch, roh_min_l, 
                  starts_pos=[], ends_pos=[]):
        """Create and returndthe hapROH dataframe."""
        
        df_full = pd.DataFrame({'Start': starts, 'End': ends,
                                'StartM': starts_map, 'EndM': ends_map, 'length': l,
                                'lengthM': l_map, 'iid': iid, "ch": ch})
        ### Add Physical Positions if given
        if len(starts_pos)>0:
            df_full["StartPosGRCh37"] = starts_pos
            df_full["EndPosGRCh37"] = ends_pos    
            
        df = df_full[df_full["lengthM"] > roh_min_l]  # Cut out long blocks
        return df

    def call_roh(self, ch=0, iid=""):
        """Call ROH of Homozygosity from Posterior Data
        bigger than cutoff
        log: Whether Posterior is given in log space"""
        r_map, pos, posterior0 = self.load_data()
        roh_post = self.modify_posterior0(posterior0)
        roh = roh_post > self.cutoff_post

        if self.output == True:
            frac_roh = np.mean(roh)
            print(f"Fraction Markers in ROH: {frac_roh:.4f}")

        # Identify Stretches by difference (up and down)
        x1 = np.hstack([[False], roh, [False]]).astype("int")  # padding
        d = np.diff(x1)
        starts = np.where(d == 1)[0]
        ends = np.where(d == -1)[0]
        l = ends - starts
        
        ### Prepare Map positions
        ends_map = r_map[ends - 1]  # -1 to stay within bounds
        starts_map = r_map[starts]
        l_map = ends_map - starts_map
        
        ### Prepare physical positions
        ends_pos = pos[ends - 1]  # -1 to stay within bounds
        starts_pos = pos[starts]
        
        # Create hapROH Dataframe
        df = self.create_df(starts, ends, starts_map, ends_map, 
                            l, l_map, iid, ch, 
                            roh_min_l=self.roh_min_l,
                            starts_pos=starts_pos, ends_pos=ends_pos)

        # Merge Blocks in Postprocessing Step
        if self.merge == True:
            df = self.merge_called_blocks(df)

        if self.snps_extend != 0:   # Extend by Nr of SNPs if needed
            df = self.snp_extend(df, r_map)

        if self.output == True:
            print(f"Called n={len(df)} ROH Blocks > {self.roh_min_l * 100} cM")
            l = np.max(df["lengthM"])
            print(f"Longest Block: {l *100:.2f} cM")

        self.df = df
        if self.save == True:
            save_folder = os.path.join(self.folder, "roh.csv")
            df.to_csv(save_folder, index=False)

            if self.output == True:
                print(f"Successfully saved to {save_folder}")

        return df

    def clean_up(self, full=True):
        """Removes all additional Data other than the
        ROH Calls and the ROH Ground Truth (To save space)"""
        keep_files = ["roh.csv", "roh_gt.csv"]
        folder = self.folder

        for the_file in os.listdir(folder):   # Walk through the Files
            file_path = os.path.join(folder, the_file)

            if os.path.isfile(file_path) and (the_file not in keep_files):
                os.unlink(file_path)   # Delete the File

                
#######################################################
class PostProcessingX(PostProcessing):
    """Class to post-process IBD on the X of 
    two males. Only difference: Two iids,
    which will get stored seperately."""
    
    def create_df(self, starts, ends, starts_map, ends_map, 
                  l, l_map, iid, ch, roh_min_l,
                  starts_pos=[], ends_pos=[]):
        """Create and returndthe hapROH dataframe.
        Difference: Here it is a IBD, so two iids are saved."""
        assert(len(iid)==2) # Sanity check
        
        df_full = pd.DataFrame({'Start': starts, 'End': ends,
                                'StartM': starts_map, 'EndM': ends_map, 'length': l,
                                'lengthM': l_map, 
                                'iid1': iid[0], "iid2":iid[1], "ch": ch})
        ### Add Physical Positions if given
        if len(starts_pos)>0:
            df_full["StartPosGRCh37"] = starts_pos
            df_full["EndPosGRCh37"] = ends_pos    
        
        
        df = df_full[df_full["lengthM"] > roh_min_l]  # Cut out long blocks
        return df
    
#######################################################


class MMR_PostProcessing(PostProcessing):
    """Class that does PostProcessing of HAPSBURG output.
    Same as PostProcessing but load Posterior differently"""

    def modify_posterior0(self, posterior0):
        """Load and return the posterior. Don't do anything"""
        roh_post = posterior0
        return roh_post

#######################################################
#######################################################


def load_Postprocessing(folder="", method="Standard", output=True, save=True):
    """Factory Method for PostProcessing class"""
    if method == "Standard":
        pp = PostProcessing(folder, output=output, save=save)
    elif method == "MMR":
        pp = MMR_PostProcessing(folder, output=output, save=save)
    elif method == "IBD_X":
        pp = PostProcessingX(folder, output=output, save=save)
    else:
        raise RuntimeError(f"Postprocessing method {method} not available!")
    return pp


#######################################################
# Do testing
if __name__ == "__main__":
    # d05e e: Error Introduced. d05: Downsampled
    folder = "./Simulated/1000G_Mosaic/TSI/ch3_10cm/output/iid0/chr3/"
    pp = PostProcessing(folder=folder)
    pp.set_params(snps_extend=0, save=False)
    df = pp.call_roh()
    print(df)

    pp.set_params(snps_extend=6, save=False)
    df = pp.call_roh()
    print(df)

"""
Class for calling ROH from Posterior Data. Saves results as a .csv
created by Pandas.
Contains Sub-Classes, as well as factory Method.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import numpy as np
import pandas as pd
import os


class PostProcessing(object):
    """Class for transition probabilities.
    Has methods to return them"""
    folder = ""         # The Folder to operate in
    r_map = []          # The Genetic Map
    posterior = []      # The Posterior
    roh_df = []        # Dataframe for Runs of Homozygosity

    cutoff = 0.8  # Cutoff Probability for ROH State
    l_cutoff = 0.01  # Cutoff in Morgan

    output = True
    save = True  # Whether to save output into Folder

    def __init__(self, folder="", load=True, output=True, save=True):
        """Initialize Class"""
        if load == True:
            self.folder = folder
            self.load_data()

        self.output = output
        self.save = save

    def load_data(self, folder=""):
        """Load and return genetic Map and Posterior"""
        if len(folder) == 0:
            folder = self.folder  # Use the Folder of the Class

        # Linkage Map
        map_path = folder + "map.csv"

        if os.path.exists(map_path):
            self.r_map = np.loadtxt(
                map_path, dtype="float", delimiter=",")
        else:
            # Eventually: Runtime Warning
            print("No Genetic Map found. Defaulting...")
            self.r_map = np.arange(np.shape(self.posterior)[1])

        # Load Posterior
        post_path = folder + "posterior.csv"
        self.posterior = np.loadtxt(post_path, dtype="float", delimiter=",")

        assert(len(self.r_map) == np.shape(self.posterior)[1])  # Sanity Check

        print(f"Successfully loaded for PP. from {folder}")

    def call_roh(self):
        """Call ROH of Homozygosity from Posterior Data"""
        roh_post = 1 - \
            np.exp(self.posterior[0, :])  # Calculate the ROH Posterior

        roh = roh_post > self.cutoff

        if self.output == True:
            frac_roh = np.mean(roh)
            print(f"Fraction Markers in ROH: {frac_roh:.4f}")

        # Identify Stretches:
        x1 = np.hstack([[False], roh, [False]]).astype("int")  # padding
        d = np.diff(x1)
        starts = np.where(d == 1)[0]
        ends = np.where(d == -1)[0]
        l = ends - starts

        ends_map = self.r_map[ends]  # Transform to Map Lengths
        starts_map = self.r_map[starts]
        l_map = ends_map - starts_map

        full_df = pd.DataFrame({'Start': starts, 'End': ends,
                                'StartM': starts_map, 'EndM': ends_map, 'length': l, 'lengthM': l_map})
        df = full_df[full_df["lengthM"] > self.l_cutoff]

        if self.output == True:
            print(f"Called n={len(df)} ROH Blocks > {self.l_cutoff * 100} cM")
            l = np.max(df["lengthM"])
            print(f"Longest Block: {l *100:.3f}")

        self.df = df
        if self.save == True:
            save_folder = self.folder + "roh.csv"
            df.to_csv(save_folder, index=False)

            if self.output==True:
                print(f"Successfully saved to {save_folder}")

        return df

#######################################################
#######################################################


def give_Postprocessing(folder="", method="Standard", output=True, save=True):
    """Factory Method for PostProcessing class"""
    pp = PostProcessing(folder, output=output, save=save)
    return pp

# Do testing
if __name__ == "__main__":
    # d05e e: Error Introduced. d05: Downsampled
    folder = "./Empirical/Sard100_0-10kROH/"
    pp = PostProcessing(folder=folder)
    pp.call_roh()

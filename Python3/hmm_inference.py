import numpy as np
import matplotlib.pyplot as plt
import os  # For Saving to Folder


#################################
#################################


class HMM_Analyze(object):
    """Analyze Class for HMMs.
    This is the main Class, all specific Inference schemes inherit from it
    and overwrite functions.
    Contains Parameters"""
    l = 1000  # Nr of the Observations
    ref_states = []  # Ref. Array of k Reference States to Copy from. [kxl]
    ob_stat = []  # The observed State [l]

    def __init__(self):
        """Initialize Class"""
        pass

    def load_data(self, folder="./Simulated/Example0/"):
        """Loads the Data"""
        self.ref_stats = np.loadtxt(
            folder + "refs.csv", dtype="int", delimiter=",")
        self.ob_stat = np.loadtxt(
            folder + "hap.csv", dtype="int", delimiter=",")

        print(f"Successfully loaded from: {folder}")


###############################
###############################
# Do some Testing
hmm = HMM_Analyze()
hmm.load_data(folder="./Simulated/Example0/")
print(np.shape(hmm.ref_stats))

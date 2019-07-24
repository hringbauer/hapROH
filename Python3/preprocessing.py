"""
Class for preprocessing and loading the Data.
Is the interface to the Folders
Has Potential for Sub-Classes for different Types of Data, as well as factory Method.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import allel  # Scikit Allel
import h5py   # For Processing HDF5s
import numpy as np
import pandas as pd
import os   # For creating folders

# Write General PreProcessing Class:
# Inherit one for real HDF5 Dataset: PreProcessingHDF5
# Inherit one for simulated Data.


class PreProcessing(object):
    """Class for PreProcessing the Data.
    Standard: Intersect Reference Data with Individual Data
    Return the Intersection Dataset
    """
    ref_folder = ""
    ind_folder = ""

    iid = ""     # Which Individual to Analyze
    ch = 0       # Which Chromosome to analyze

    def __init__(self):
        """Initialize Class"""
        raise NotImplementedError()

    def load_data(self, iid="MA89", ch=6, n_ref=503):
        """Return Refererence Matrix [k,l], Genotype/Readcount Matrix [2,l]
        as well as linkage map [l] """
        raise NotImplementedError()

    def set_params(self, **kwargs):
        """Set the Parameters.
        Takes keyworded arguments"""
        for key, value in kwargs.items():
            setattr(self, key, value)


class PreProcessingHDF5(PreProcessing):
    """Class for PreProcessing the Data.
    Standard: Intersect Reference Data with Individual Data
    Return the Intersection Dataset
    """
    out_folder = ""    # Where to Save to (Returned to main program)
    h5_path_targets = "./../ancient-sardinia/output/h5_rev/mod_reich_sardinia_ancients_rev_mrg_dedup_3trm_anno.h5"
    meta_path_targets = "./../ancient-sardinia/output/meta/meta_rev_final.csv"

    # Path of 1000G (without chromosome part)
    h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/Eur1240chr"
    meta_path_ref = "./Data/1000Genomes/Individuals/meta_df.csv"
    excluded = ["TSI", ]  # List of excluded Populations in Meta

    prefix_out_data = ""  # Prefix of the Outdata (should be of form "path/")

    save = True
    output = True
    readcounts = False   # Whether to return Readcounts
    diploid_ref = True   # Whether to use diploid Reference Individuals
    destroy_phase = True  # Whether to destroy Phase of Target Individual
    only_calls = True    # Whether to downsample to markers with no missing data

    def __init__(self, save=True, output=True):
        """Initialize Class.
        Ind_Folder: Where to find individual
        iid & chr: Individual and Chromosome.
        save: """
        self.save = save
        self.output = output
        # print(os.getcwd()) # Show the current working directory

    def set_output_folder(self, iid, ch):
        """Set the output folder."""
        out_folder = "./Empirical/1240k/" + \
            str(iid) + "/chr" + str(ch) + "/" + self.prefix_out_data
        return out_folder

    def get_index_iid(self, iid, fs=0):
        """Get the Index of IID in fs
        iid to extract. fs reference HDF5"""
        meta_df = pd.read_csv(self.meta_path_targets)
        assert(len(meta_df) == np.shape(fs["calldata/GT"])[1])  # Sanity Check

        id_obs = np.where(meta_df["iid"] == iid)[0]
        if len(id_obs) == 0:
            raise RuntimeError(f"Individual {iid} not found in {self.meta_path_targets}!")
        return id_obs[0]

    def get_ref_ids(self, f, n_ref):
        """OVERWRITE: Get the Indices of the individuals
        in the HDF5 to extract. Here: Allow to subset for Individuals from
        different 100G Populations"""

        # Load Meta Population File
        meta_df = pd.read_csv(self.meta_path_ref, sep="\t")

        iids = np.where(~meta_df["pop"].isin(self.excluded))[0]
        print(f"{len(iids)} / {len(meta_df)} Individuals included in Reference")
        return iids[:n_ref]   # Return # n_ref Individual Indices

    def load_data(self, iid="MA89", ch=6, n_ref=503, folder=""):
        """Return Matrix of reference [k,l], Matrix of Individual Data [2,l],
        as well as linkage Map [l]"""

        if self.output==True:
            print(f"Loading Individual: {iid}")

        # Attach Part for the right Chromosome
        h5_path1000g = self.h5_path1000g + str(ch) + ".hdf5"

        # Def Set the output folder:
        out_folder = self.set_output_folder(iid, ch)

        # Create Output Folder if needed
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)

        # Set important "Steady Paths":
        h5_path_targets = self.h5_path_targets

        # Load and Merge the Data
        fs = self.load_h5(h5_path_targets)
        f1000 = self.load_h5(h5_path1000g)
        i1, i2 = self.merge_2hdf(fs, f1000)

        id_obs = self.get_index_iid(iid, fs)
        ids_ref = self.get_ref_ids(f1000, n_ref)

        # All 503 EUR Samples as Reference (first Chromosome)
        markers = np.arange(0, len(i1))  # Which Markers to Slice out
        markers_obs = i1[markers]
        markers_ref = i2[markers]

        gts_ind, gts, read_counts, r_map = self.extract_snps(out_folder, f1000, fs, ids_ref, id_obs,
                                                             markers_ref, markers_obs)

        if self.only_calls == True:
            called = self.markers_called(gts_ind, read_counts)
            gts_ind = gts_ind[:, called]
            gts = gts[:, called]
            r_map = r_map[called]
            read_counts = read_counts[:, called]
            if self.output == True:
                print(f"Markers called {np.sum(called)} / {len(called)}")

        if self.save == True:
            self.save_info(out_folder, r_map,
                           gt_individual=gts_ind, read_counts=read_counts)

        if self.readcounts == True:   # Switch to Readcount
            if self.output == True:
                print(f"Loading Readcounts...")
                print(f"Mean Readcount markers loaded: {np.mean(read_counts) * 2:.5f}")
            gts_ind = read_counts

        if self.destroy_phase == True:     # Destroy Phase
            if self.output == True:
                print("Shuffling phase of target...")
            gts_ind = self.destroy_phase(gts_ind)

        return gts_ind, gts, r_map, out_folder

     ################################################
     # Some Helper Functions

    def destroy_phase(self, gts_ind):
        """Randomly shuffles phase for gts [2,n_loci]"""
        assert(np.shape(gts_inds)[0] == 2)

        n_loci = np.shape(gts_ind)[1]
        phases = np.random.randint(2, size=n_loci)  # Do the random shuffling

        gts_ind_new = np.zeros(np.shape(gts_inds), dtype="int")
        gts_ind_new[0, :] = gts_ind[phases, np.arange(n_loci)]
        gts_ind_new[1, :] = gts_ind[1 - phases, np.arange(n_loci)]
        return gts_ind_new

    def markers_called(self, gts_ind, read_counts):
        """Return boolean array of markers which are called"""
        called = []
        if self.readcounts == False:
            called = (gts_ind[0, :] > -1)  # Only Markers with calls
        elif self.readcounts == True:
            read_depth = np.sum(read_counts, axis=0)
            called = read_depth > 0
        else:
            raise RuntimeError("Invalid Mode")
        return called

    def load_h5(self, path):
        """Load and return the HDF5 File from Path"""
        f = h5py.File(path, "r")  # Load for Sanity Check. See below!
        print("\nLoaded %i variants" % np.shape(f["calldata/GT"])[0])
        print("Loaded %i individuals" % np.shape(f["calldata/GT"])[1])
        # print(list(f["calldata"].keys()))
        # print(list(f["variants"].keys()))
        print(f"HDF5 loaded from {path}")
        return f

    def merge_2hdf(self, f, g, ch=1):
        """ Merge two HDF 5 f and g. Return Indices of Overlap Individuals.
        f is Sardinian HDF5,
        g the Reference HDF5
        ch: Integer, which Chromosome to use"""

        pos1 = f["variants/POS"]
        pos2 = g["variants/POS"]

        # Check if in both Datasets
        b, i1, i2 = np.intersect1d(pos1, pos2, return_indices=True)

        # Sanity Check if Reference is the same
        ref1 = np.array(f["variants/REF"])[i1]
        ref2 = np.array(g["variants/REF"])[i2]
        alt1 = np.array(f["variants/ALT"])[i1]
        alt2 = np.array(g["variants/ALT"])[i2, 0]

        # Downsample to Site where both Ref and Alt are the same
        same = (ref1 == ref2)

        both_same = (ref1 == ref2) & (alt1 == alt2)
        i11 = i1[both_same]
        i22 = i2[both_same]

        if self.output == True:
            print(f"\nIntersection on Positions: {len(b)}")
            print(f"Nr of Matching Refs: {np.sum(same)} / {len(same)}")
            print(f"Full Intersection Ref/Alt Identical: {len(i11)} / {len(both_same)}")

        return i11, i22

    def save_info(self, folder, cm_map, gt_individual=[], read_counts=[]):
        """Save Linkage Map, Readcount and Genotype Data per Individual.
        (Needed for latter Plotting)
        Genotypes Individual: If given, save as well"""

        # Save the cmap
        np.savetxt(folder + "map.csv", cm_map, delimiter=",",  fmt='%.8f')

        if len(gt_individual) > 0:
            np.savetxt(folder + "hap.csv", gt_individual,
                       delimiter=",",  fmt='%i')
        if len(read_counts) > 0:
            np.savetxt(folder + "readcounts.csv", read_counts,
                       delimiter=",",  fmt='%i')

        if self.output == True:
            print(f"Successfully saved to: {folder}")

    #######################################
    # Code for saving Haplotype
    def extract_snps(self, folder, ref_hdf5, obs_hdf5, ids_ref, id_obs,
                     marker_ref, marker_obs):
        """Save Folder with all relevant Information.
        Folder: Where to save to
        ref_hdf5: Reference HDF5
        obs_hdf5: Observed HDF5
        ids_ref: Indices of reference Individuals to save
        ids_obs: Indices of observed Individuals
        marker_ref: Indices of reference Markers
        marker_obs: Indices of observed Markers
        only_calls: Whether to Only Include Markers with Calls"""
        assert(len(marker_ref) == len(marker_obs)
               )  # If reference and observe dataset are the same

        # Extract Reference Individuals (first haplo)
        gts = ref_hdf5["calldata/GT"][:, ids_ref, 0]  # Only first IID
        gts = gts[marker_ref, :].T       # Important: Swap of Dimensions!!

        if self.diploid_ref == True:   # In case diploid reference Samples
            gts1 = ref_hdf5["calldata/GT"][:, ids_ref, 1]  # The second allele
            # Important: Swap of Dimensions!!
            gts1 = gts1[marker_ref, :].T
            gts = np.concatenate((gts, gts1), axis=0)  # Add two dataframes

        if self.output == True:
            print(f"Extraction of {len(gts)} Haplotypes Complete!")

        # Extract target individual Genotypes
        gts_ind = obs_hdf5["calldata/GT"][:, id_obs, :]
        gts_ind = gts_ind[marker_obs, :].T

        # Extract Readcounts
        read_counts = obs_hdf5["calldata/AD"][:, id_obs, :]
        read_counts = read_counts[marker_obs, :].T

        # Extract Linkage map
        r_map = np.array(ref_hdf5["variants/MAP"]
                         )[marker_ref]  # Load the LD Map

        # np.savetxt(folder + "ind.csv", [id_obs], delimiter=",",  fmt='%i')
        # Return Genotypes/Readcounts Individual, Genotypes Reference and Recombination Map
        return gts_ind, gts, read_counts, r_map


###########################################

class PreProcessingHDF5Sim(PreProcessingHDF5):
    """Class for PreProcessing simulated 1000 Genome Data (Mosaic Simulations).
    Same as PreProcessingHDF5 but with the right Folder Structure
    MODIFY
    """

    out_folder = ""    # Where to save to
    prefix_out_data = ""  # Prefix of the Outdata
    meta_path_targets = "./../ancient-sardinia/output/meta/meta_final.csv"
    h5_folder = ""    # The H5 Folder
    h5_path_targets = ""
    # Path of 1000G (without chromosome part):
    h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/Eur1240chr"

    def set_output_folder(self, iid, ch):
        """Set the output folder of where to save the result to."""

        out_folder = self.h5_folder + "output/" + \
            str(iid) + "/chr" + str(ch) + "/" + self.prefix_out_data

        return out_folder

    def get_index_iid(self, iid, fs=0):
        """OVERWRITE: Get the Index of IID in fs
        iid to extract. fs reference HDF5.
        Difference: Here found directly in HDF Samples"""

        iids = np.array(fs['samples'])
        assert(len(iids) == np.shape(fs["calldata/GT"])[1])  # Sanity Check

        id_obs = np.where(iids == iid)[0]

        if len(id_obs) == 0:
            raise RuntimeError(f"Individual {iid} not found in {self.meta_path_targets}!")
        return id_obs[0]

    def set_folder(self, folder_path):
        """Method to manually set the Input and Output Folder"""
        self.h5_folder = folder_path
        self.h5_path_targets = folder_path + "data.h5"

    def set_prefix_out_data(self, prefix):
        """Modify the Prefix of the Output-File
        LEGACY: Should be done via set_params()"""
        self.prefix_out_data = prefix

    def set_exclude_pops(self, pops=["TSI", ]):
        """Method to manually set the excluded populations.
        LEGACY: Should be done via set_params()"""
        self.excluded = pops

############################################
############################################


class PreProcessingFolder(PreProcessing):
    """Preprocessing if data has been saved into a folder
    (such as in a Simulation)
    """
    save = True
    output = True
    readcounts = False   # Whether to return Readcounts

    def __init__(self, save=True, output=True):
        """Initialize Class.
        Ind_Folder: Where to find individual
        iid & chr: Individual and Chromosome.
        save: """
        self.save = save
        self.output = output

    def load_data(self, iid="MA89", ch=6, n_ref=503, folder=""):
        """Return Matrix of reference [k,l], Matrix of Individual Data [2,l],
        as well as linkage Map [l]"""

        gts = np.loadtxt(
            folder + "refs.csv", dtype="int", delimiter=",")

        gts_ind = np.loadtxt(
            folder + "hap.csv", dtype="int", delimiter=",")

        n_snps = np.shape(gts_ind)[1]
        r_map = self.load_linkage_map(folder, n_snps)

        if self.output:
            print(f"Successfully loaded Data from: {folder}")

        return gts_ind, gts, r_map, folder

    def load_linkage_map(self, folder, nr_snps):
        """Load and Return the Linkage Map"""
        map_path = folder + "map.csv"

        if os.path.exists(map_path):
            r_map = np.loadtxt(
                map_path, dtype="float", delimiter=",")

        else:
            # Eventually: Runtime Warning
            print("No Genetic Map found. Defaulting...")
            r_map = np.arange(nr_snps)

        assert(len(r_map) == nr_snps)  # Sanity Check

        return r_map


############################################
############################################


############################################
############################################
# Do a Factory Method that can be imported.


def load_preprocessing(p_model="SardHDF5", save=True, output=True):
    """Load the Transition Model"""

    if p_model == "SardHDF5":
        p_obj = PreProcessingHDF5(save=save, output=output)

    elif p_model == "MosaicHDF5":
        p_obj = PreProcessingHDF5Sim(save=save, output=output)

    elif p_model == "Folder":
        p_obj = PreProcessingFolder(save=save, output=output)

    else:
        raise NotImplementedError("Transition Model not found!")

    return p_obj


# For testing the Module
if __name__ == "__main__":
    pp = load_preprocessing(p_model="Folder", save=False, output=True)
    gts_ind, gts, r_map, out_folder = pp.load_data(
        iid="MA89", ch=3, n_ref=503, folder="./Simulated/Test20r/")
    print(gts_ind[:2, :4])
    print(np.shape(gts_ind))
    print(r_map[:5])
    print(np.shape(r_map))
    print(gts[:10, :2])
    print(np.shape(gts))
    print(out_folder)

"""
Class for preprocessing and loading the Data.
Is the interface to Data Folders
Has Sub-Classes for different Types of Data, as well as factory Method
that returns the right subclass based on a keyword. Use this 
factory method to load classes.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

#import allel  # Scikit Allel
import h5py   # For Processing HDF5s
import numpy as np
import pandas as pd
import os   # For creating folders


#Assume hapsburg directory is in root
from hapsburg.PackagesSupport.loadEigenstrat.loadEigenstrat import load_eigenstrat

# Write General PreProcessing Class: (PreProcessing)
# Inherit one for real HDF5 Dataset: PreProcessingHDF5
# > Inherit one for simulated Data: PreProcessingHDF5Sim
# > Inherit one for Eigenstrat Data: PreProcessingEigenstrat


class PreProcessing(object):
    """Class for PreProcessing the Data.
    Standard: Intersect Reference Data with Individual Data
    Return the Intersection Dataset
    """
    save = True
    output = True

    iid = ""     # Which Individual to Analyze
    ch = 0       # Which Chromosome to analyze
    segM = []    # Segment of Chromosome to analyze (in Morgan)
    n_ref = 0   # Nr of diploid reference Individuals to use (diploid)

    def __init__(self):
        """Initialize Class"""
        raise NotImplementedError()

    def load_data(self, iid="MA89", ch=6):
        """Return Refererence Matrix [k,l], Genotype/Readcount Matrix [2,l]
        as well as linkage map [l] """
        raise NotImplementedError()

    def set_params(self, **kwargs):
        """Set the Parameters.
        Takes keyworded arguments"""
        for key, value in kwargs.items():
            setattr(self, key, value)

    def set_output_folder(self, iid, ch):
        """Set the output folder after folder_out.
        General Structure for HAPSBURG: folder_out/iid/chrX/"""
        out_folder = os.path.join(self.folder_out, str(
            iid),  "chr" + str(ch), self.prefix_out_data)

        if not os.path.exists(out_folder):   # Create Output Folder if needed
            if self.output == True:
                print(f"Creating folder {out_folder}...")
            os.makedirs(out_folder)

        return out_folder

class PreProcessingHDF5(PreProcessing):
    """Class for PreProcessing the Data.
    Standard: Intersect Reference Data with Individual Data
    Return the Intersection Dataset
    """
    path_targets = "./../ancient-sardinia/output/h5_rev/mod_reich_sardinia_ancients_rev_mrg_dedup_3trm_anno.h5"
    #meta_path_targets = "./../ancient-sardinia/output/meta/meta_rev_final.csv" Meta now encoded in f["samples"]!!

    # Path of 1000G (without chromosome part)
    h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/Eur1240chr"
    meta_path_ref = "./Data/1000Genomes/Individuals/meta_df.csv"
    excluded = ["TSI", ]  # List of excluded Populations in Meta

    folder_out = "./Empirical/1240k/"  # Base Path of the Output Folder
    prefix_out_data = ""  # Prefix of the Outdata (should be of form "path/")

    readcounts = False    # Whether to return Readcounts
    diploid_ref = True    # Whether to use diploid Reference Individuals
    random_allele = True  # Whether to pick one of two alleles at Target Individual at random
    only_calls = True     # Whether to downsample to markers with no missing data
    max_mm_rate = 0.9     # Maximal mismatch rate ref/alt alleles between target and ref

    def __init__(self, save=True, output=True):
        """Initialize Class.
        Ind_Folder: Where to find individual
        iid & chr: Individual and Chromosome.
        save: """
        self.save = save
        self.output = output
        # print(os.getcwd()) # Show the current working directory

    def get_index_iid_legacy(self, iid, fs=0):
        """Get the Index of IID in fs
        iid to extract. fs reference HDF5"""
        meta_df = pd.read_csv(self.meta_path_targets)
        assert(len(meta_df) == np.shape(fs["calldata/GT"])[1])  # Sanity Check

        id_obs = np.where(meta_df["iid"] == iid)[0]
        if len(id_obs) == 0:
            raise RuntimeError(f"Individual {iid} not found in {self.meta_path_targets}!")
        return id_obs[0]
    
    def get_index_iid(self, iid, f=0, samples_field="samples"):
        """Get the Index of IID in fs
        iid to extract. fs reference HDF5"""
        assert(len(f[samples_field]) == np.shape(f["calldata/GT"])[1])  # Sanity Check

        id_obs = np.where(f[samples_field][:] == iid)[0]
        if len(id_obs) == 0:
            raise RuntimeError(f"Individual {iid} not found in H5 field {samples_field}!")
        return id_obs[0]

    def get_ref_ids(self, f):
        """OVERWRITE: Get the Indices of the individuals
        in the HDF5 to extract. Here: Allow to subset for Individuals from
        different 100G Populations"""

        # Load Meta Population File
        meta_df = pd.read_csv(self.meta_path_ref, sep="\t")
        assert(len(meta_df) == np.shape(f["calldata/GT"])[1])  # Sanity Check

        iids = np.where(~meta_df["pop"].isin(self.excluded))[0]
        if self.output:
            print(f"{len(iids)} / {len(meta_df)} Individuals included in Reference")
            print(f"Extracting up to {self.n_ref} Individuals")
        return iids[:self.n_ref]   # Return up to n_ref Individual Indices

    def load_data(self, iid="MA89", ch=6):
        """Return Matrix of reference [k,l], Matrix of Individual Data [2,l],
        as well as linkage Map [l]"""

        if self.output == True:
            print(f"Loading Individual: {iid}")

        # Attach Part for the right Chromosome
        h5_path1000g = self.h5_path1000g + str(ch) + ".hdf5"

        # Def Set the output folder:
        out_folder = self.set_output_folder(iid, ch)

        # Load and Merge the Data
        fs = self.load_h5(self.path_targets)
        f1000 = self.load_h5(h5_path1000g)
        i1, i2 = self.merge_2hdf(fs, f1000)

        id_obs = self.get_index_iid(iid, fs)
        ids_ref = self.get_ref_ids(f1000)

        # All 503 EUR Samples as Reference (first Chromosome)
        markers = np.arange(0, len(i1))  # Which Markers to Slice out
        markers_obs = i1[markers]
        markers_ref = i2[markers]

        # Load Target Dataset
        gts_ind = self.extract_snps_hdf5(
            fs, [id_obs], markers_obs, diploid=True)
        read_counts = self.extract_rc_hdf5(fs, id_obs, markers_obs)

        # Load Reference Dataset
        gts = self.extract_snps_hdf5(
            f1000, ids_ref, markers_ref, diploid=self.diploid_ref)
        r_map = self.extract_rmap_hdf5(f1000, markers_ref)  # Extract LD Map

        # Do optional Processing Steps (based on boolean flags in class)
        gts_ind, gts, r_map, out_folder = self.optional_postprocessing(
            gts_ind, gts, r_map, out_folder, read_counts)

        return gts_ind, gts, r_map, out_folder

    def optional_postprocessing(self, gts_ind, gts, r_map, out_folder, read_counts=[]):
        """Postprocessing steps of gts_ind, gts, r_map, and the folder,
        based on boolean fields of the class."""

        if self.only_calls == True:
            called = self.markers_called(gts_ind, read_counts)
            gts_ind = gts_ind[:, called]
            gts = gts[:, called]
            r_map = r_map[called]
            if len(read_counts) > 0:
                read_counts = read_counts[:, called]
            if self.output:
                print(f"Reduced to markers called {np.shape(gts)[1]} / {len(called)}")
                print(f"Fraction SNPs covered: {np.shape(gts)[1] / len(called):.4f}")

        if self.save == True:
            self.save_info(out_folder, r_map,
                           gt_individual=gts_ind, read_counts=read_counts)

        if (self.readcounts == True) and len(read_counts) > 0:   # Switch to Readcount
            if self.output:
                print(f"Loading Readcounts...")
                print(f"Mean Readcount markers loaded: {np.mean(read_counts) * 2:.5f}")
            gts_ind = read_counts
        
        ### Shuffle Target Allele     
        if self.random_allele == True:     
            if self.output == True:
                print("Shuffling phase of target...")
            gts_ind = self.destroy_phase_func(gts_ind)

        return gts_ind, gts, r_map, out_folder

     ################################################
     # Some Helper Functions

    def destroy_phase_func(self, gts_ind):
        """Randomly shuffles phase for gts [2,n_loci]"""
        assert(np.shape(gts_ind)[0] == 2)

        n_loci = np.shape(gts_ind)[1]
        phases = np.random.randint(2, size=n_loci)  # Do the random shuffling

        gts_ind_new = np.zeros(np.shape(gts_ind), dtype="int")
        gts_ind_new[0, :] = gts_ind[phases, np.arange(n_loci)]
        gts_ind_new[1, :] = gts_ind[1 - phases, np.arange(n_loci)]
        return gts_ind_new

    def markers_called(self, gts_ind, read_counts):
        """Return boolean array of markers which are called"""
        called = []
        if (self.readcounts == False) or (len(read_counts) == 0):
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
        
    def merge_2hdf(self, f, g):
        """ Merge two HDF 5 f and g. Return Indices of Overlap Individuals.
        f is Sardinian HDF5,
        g the Reference HDF5"""

        pos1 = f["variants/POS"]
        pos2 = g["variants/POS"]

        # Check if in both Datasets
        b, i1, i2 = np.intersect1d(pos1, pos2, return_indices=True)

        ### Sanity Check if Reference is the same
        ref1 = np.array(f["variants/REF"])[i1]
        ref2 = np.array(g["variants/REF"])[i2]
        alt1 = np.array(f["variants/ALT"])[i1]
        
        if len(np.shape(g["variants/ALT"]))>1:   # If multiple ALT Alleles
            alt2 = np.array(g["variants/ALT"])[i2, 0]
        else:
            alt2 = np.array(g["variants/ALT"])[i2]

        ### Downsample to Site where both Ref and Alt are the same
        same = (ref1 == ref2)
        both_same = (ref1 == ref2) & (alt1 == alt2)
        i11 = i1[both_same]
        i22 = i2[both_same]

        if self.output == True:
            print(f"Nr of Matching Refs: {np.sum(same)} / {len(same)} SNPs")
            print(f"Both Ref/Alt Identical: {len(i11)} / {len(both_same)}")
        
        mm_frac = (len(i11)/len(both_same))
        if mm_frac <= self.max_mm_rate:
            error = f"Ref/Alt Base Match fraction {mm_frac:.4f} < {self.max_mm_rate}"
            error1 = " Check genetic map" # Maybe enter what desired version should be
            raise RuntimeError(error + error1)

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

    def extract_snps_hdf5(self, h5, ids, markers, diploid=False, dtype=np.int8):
        """Extract from h5 on ids and markers.
        If diploid, concatenate haplotypes"""
        # Extract Reference Individuals (first haplo)
        gts = h5["calldata/GT"][:, ids, 0].astype(dtype)  # Only first IID
        gts = gts[markers, :].T       # Important: Swap of Dimensions!!

        if diploid == True:   # In case diploid: Second Alllele
            gts1 = h5["calldata/GT"][:, ids, 1].astype(dtype)
            gts1 = gts1[markers, :].T
            gts = np.concatenate((gts, gts1), axis=0)  # Combine dataframes

        if self.output == True:
            print(f"Extraction of {len(gts)} Haplotypes complete")
        return gts

    def extract_rc_hdf5(self, h5, ids, markers, dtype=np.int16):
        """Extract Readcount data from from h5 on single (!) id and markers"""
        read_counts = h5["calldata/AD"][:, ids, :].astype(dtype)
        read_counts = read_counts[markers, :].T
        return read_counts

    def extract_rmap_hdf5(self, h5, markers):
        """Get the Linkage Map from h5"""
        r_map = np.array(h5["variants/MAP"])[markers]
        return r_map

    def get_segment(self, r_map, markers_obs, markers_ref):
        """Extract only markers in self.segM and downsamples
        the key array. Return downsampled versions."""
        segM = self.segM
        assert(len(seg) == 2)  # Sanity Check
        in_seg = (r_map > seg[0]) & (r_map < segM[1])
        if self.output == True:
            print(f"Extracting {np.sum(in_seg)}/{len(in_seg)} SNPs"
                   "within {segM[0]}-{segM[1]}")
        return r_map[in_seg], markers_obs[in_seg], markers_ref[in_seg]

###########################################


class PreProcessingEigenstrat(PreProcessingHDF5):
    """Class for PreProcessing Eigenstrat Files
    Same as PreProcessingHDF5 for reference, but with Eigenstrat coe
    for target
    """
    #meta_path_targets = "" # Meta not needed any more, only .ind file
    path_targets = ""  # Base Path of the Eigenstrat
    # Path of 1000G (without chromosome part):
    h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/Eur1240chr"
    packed = True  # Whether to use packed or unpacked Eigenstrat
    sep = r"\s+"   # Which Column Separator to use in ind and snp File
    flipstrand = True # Flip Strand if both alleles matching, but flipped
    
    def __init__(self, save=True, output=True, packed=True, sep= r"\s+"):
        """Initialize Class.
        save: Whether to save the data
        output: Whether to print output
        packed: Whether Eigenstrat .geno file is packed (binary)
        sep: What seperator to use in .ind and .geno File"""
        self.save = save
        self.output = output
        self.packed = packed
        self.sep = sep

    def load_data(self, iid="MA89", ch=6):
        """Return Matrix of reference [k,l], Matrix of Individual Data [2,l],
        as well as linkage Map [l] and the output folder.
        Save the loaded data if self.save==True
        Various modifiers in class fields (check also PreProcessingHDF5)"""
        self.ch = ch  # To remember Chromosome
        if self.output == True:
            print(f"Loading Individual: {iid}")

        out_folder = self.set_output_folder(iid, ch)

        ### Load Reference HDF5 object
        h5_path1000g = self.h5_path1000g + str(ch) + ".hdf5"
        f1000 = self.load_h5(h5_path1000g)

        ### Load Eigenstrat object
        es = load_eigenstrat(base_path=self.path_targets, 
                             sep=self.sep, packed=self.packed)

        markers_obs, markers_ref, flipped = self.merge_es_hdf5(es, f1000)

        id_obs = es.get_index_iid(iid)  # Get Index from Eigenstrat
        ids_ref = self.get_ref_ids(f1000)

        r_map = self.extract_rmap_hdf5(f1000, markers_ref)  # Extract LD Map

        if len(self.segM) > 0:
            r_map, markers_obs, markers_ref = self.get_segment(
                r_map, markers_obs, markers_ref)

        ### Do Extraction for Reference
        gts = self.extract_snps_hdf5(
            f1000, ids_ref, markers_ref, diploid=self.diploid_ref)

        ### Extraction for target (no RC here)
        gts_ind = self.extract_snps_es(es, id_obs, markers_obs)
        
        ### Flip target genotypes where flipped
        if self.flipstrand==True:
            if self.output:
                print(f"Flipping Ref/Alt in target for {np.sum(flipped)} SNPs...")
            flip_idcs = flipped & (gts_ind[0,:]>=0) # Where Flip AND Genotype Data
            gts_ind[:,flip_idcs] = 1 - gts_ind[:,flip_idcs]

        ### Do optional Processing Steps (only covered, destroy phase, save)
        gts_ind, gts, r_map, out_folder = self.optional_postprocessing(
            gts_ind, gts, r_map, out_folder)

        return gts_ind, gts, r_map, out_folder

    def extract_snps_es(self, es, id, markers):
        """Use Eigenstrat object. Extract genotypes for individual index i
        (integer) for
        list of markers. Do conversion from Eigenstrat GT to format
        used here"""
        gts_ind = es.extract_snps(id, markers, conversion=True)
        return gts_ind

    def merge_es_hdf5(self, es, f_ref):
        """Merge Eigenstrat and HDF5 Loci, return intersection indices
        es: LoadEigenstrat Object
        f_ref: Reference HDF5
        if self.flipstrand return indices [its] [its] 
        and flip boolean vector [its]"""
        pos1, idcs = es.give_positions(ch = self.ch)
        pos2 = f_ref["variants/POS"]

        # Check if in both Datasets
        b, i1, i2 = np.intersect1d(pos1, pos2, return_indices=True)

        # Load Ref/Alt alleles at intersecting Loci
        ref1, alt1 = es.give_ref_alt(ch=self.ch)
        ref1, alt1 = ref1[i1], alt1[i1]  # Subset to intersection

        ref2 = np.array(f_ref["variants/REF"])[i2]
        
        if len(np.shape(f_ref["variants/ALT"]))>1:  # If multiple alt alleles
            alt2 = np.array(f_ref["variants/ALT"])[i2, 0]
        else:
            alt2 = np.array(f_ref["variants/ALT"])[i2]

        # Downsample to Site where both Ref and Alt are identical        
        same = (ref1 == ref2)
        both_same = (ref1 == ref2) & (alt1 == alt2)
        flipped = (ref1 == alt2) & (alt1 == ref2)

        if self.output:
            print(f"\nIntersection on Positions: {len(b)}")
            print(f"Nr of Matching Refs: {np.sum(same)} / {len(same)}")
            print(f"Ref/Alt Matching: {np.sum(both_same)} / {len(both_same)}")
            print(f"Flipped Ref/Alt Matching: {np.sum(flipped)}")
            print(f"Together: {np.sum(flipped|both_same)} / {len(both_same)}")
            
        if self.flipstrand:
            i11 = idcs[i1[both_same | flipped]]  # Give the Indices in full array (all ch)
            i22 = i2[both_same | flipped]
            ### Identify Flipped SNPs indices in OR indices
            # Assume both_same and flipped have NO intersection!
            # Could clip to 1 then intersection possible
            flipped = (np.cumsum(both_same + flipped))[flipped]-1
            
        else:
            i11 = idcs[i1[both_same]]  # Give the Indices in full array (all ch)
            i22 = i2[both_same]
            flipped = [] # Numpy Indices to flip nothing
            
        flip = np.zeros(len(i11), dtype="bool")
        flip[flipped]=True
            
        mm_frac = (len(i11)/len(both_same))
        if mm_frac <= self.max_mm_rate:
            error = f"Ref/Alt Flip+Noflip fraction {mm_frac:.4f} < {self.max_mm_rate}"
            error1 = " Check genetic map" # Maybe enter what desired version should be
            raise RuntimeError(error + error1)

        return i11, i22, flip


###########################################

class PreProcessingHDF5Sim(PreProcessingHDF5):
    """Class for PreProcessing simulated 1000 Genome Data (Mosaic Simulations).
    Same as PreProcessingHDF5 but with the right Folder Structure
    MODIFY
    """

    out_folder = ""    # Where to save to
    prefix_out_data = ""  # Prefix of the Outdata
    #meta_path_targets = ""  # Not needed any more
    path_targets = ""
    # Path of 1000G (without chromosome part):
    h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/Eur1240chr"

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
# Do a Factory Method that can be imported.


def load_preprocessing(p_model="SardHDF5", save=True, output=True):
    """Load the Transition Model"""

    if p_model == "SardHDF5":
        p_obj = PreProcessingHDF5(save=save, output=output)
    elif p_model == "MosaicHDF5":
        p_obj = PreProcessingHDF5Sim(save=save, output=output)
    elif p_model == "Folder":
        p_obj = PreProcessingFolder(save=save, output=output)
    elif p_model == "EigenstratPacked":
        p_obj = PreProcessingEigenstrat(save=save, output=output,
                                        packed=True, sep=r"\s+")
    elif p_model == "EigenstratUnpacked":
        p_obj = PreProcessingEigenstrat(save=save, output=output,
                                        packed=False, sep=r"\s+")
    else:
        raise NotImplementedError(f"Preprocessing Model string {p_model} not found.")

    return p_obj


############################################
# For testing the Module
if __name__ == "__main__":
    # Test Loading Eigenstrat
    # pp = load_preprocessing(p_model="Eigenstrat", save=False, output=True)
    # pp.set_params(path_targets="./Data/ReichLabEigenstrat/Olalde2019/Olalde_et_al_genotypes",
    #            folder_out="./Empirical/ES_Test/", only_calls=True)

    # Test Loading HDF5
    pp = load_preprocessing(p_model="MosaicHDF5", save=False, output=True)
    pp.set_params(h5_path_targets="./Simulated/1000G_Mosaic/TSI5/ch3_4cm/data.h5",
                  folder_out="./Empirical/ES_Test/", only_calls=True,
                  excluded=["TSI", ])

    gts_ind, gts, r_map, out_folder = pp.load_data(
        iid="iid5", ch=3, n_ref=100)

    # gts_ind, gts, r_map, out_folder = pp.load_data(
    #    iid="MA89", ch=3, n_ref=503, folder="./Simulated/Test20r/")
    moi = range(1050, 1070)
    print(r_map[moi])
    print(np.shape(r_map))
    print(gts_ind[:2, moi])
    print(np.shape(gts_ind))
    print(gts[:5, moi])
    print(np.shape(gts))
    print(out_folder)

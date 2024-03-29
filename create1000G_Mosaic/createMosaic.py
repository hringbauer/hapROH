"""
Class for preparing 1000 Genome Mosaic Data
Loads 1000 Genome data. Creates Mosaic with copyied in ROH
Is a class that can load and save the data.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import h5py   # For Processing HDF5s
import numpy as np
import pandas as pd
import os as os
import sys


class Mosaic_1000G(object):
    """Class for Preparing 1000G Mosaics.
    Given some Parameters it creates the Mosaic individuals and returns their genotypes"""

    ch = 0  # Which Chromosome to analyze
    e_rate_ref = 1e-3
    path1000G = ""  # Path to the 1000 Genome Data
    pop_path = ""  # Path to Population Information
    save_path = ""  # Where to save the new HDF5 to

    output = True  # whether to Print Output
    f = 0          # HDF5 with 1000 Genome Data
    meta_df = 0    # Meta_df with Meta Information of the 1000 Genome Data

    def __init__(self, ch=3, e_rate_ref=1e-3, path1000G="./Data/1000Genomes/HDF5/1240kHDF5/Eur1240chr",
                 pop_path="./Data/1000Genomes/integrated_call_samples_v3.20130502.ALL.panel",
                 save_path=""):
        """ch: Which chromosome to loadself.
        pop_path: Where to load from
        save_path: Where to save the new HDF5 to"""
        print("\nStarted Mosaic Object. Working Directory:")
        print(os.getcwd()) # Show the current working directory)

        self.e_rate_ref = e_rate_ref
        # Set Path of 1000G (without chromosome part)
        self.path1000G = path1000G + str(ch) + ".hdf5"
        print(f'constructor of Mosaic_1000G: {path1000G}')
        self.pop_path = pop_path
        self.ch = ch
        # Load some Data:
        self.f = self.load_h5()
        # Load the Individual Names and merge in Population Data
        self.meta_df = self.load_meta_df()

    def load_h5(self, path=""):
        """Load and return the HDF5 File from Path"""
        if len(path) == 0:
            path = self.path1000G

        f = h5py.File(path, "r")  # Load for Sanity Check. See below!

        if self.output == True:
            print("\nLoaded %i variants" % np.shape(f["calldata/GT"])[0])
            print("Loaded %i individuals" % np.shape(f["calldata/GT"])[1])
            # print(list(f["calldata"].keys()))
            # print(list(f["variants"].keys()))
            print(f"HDF5 loaded from {path}")
            print(np.shape(f["calldata/GT"]))
        return f

    def load_meta_df(self, path=""):
        """Load, Postprocess and return the Metafile.
        Requires that f has been loaded"""
        if len(path) == 0:
            path = self.pop_path
        meta_df = pd.read_csv(path, sep="\t")

        iids = np.array(self.f["samples"]).astype("str")
        iids0 = np.char.split(iids, sep='_', maxsplit=1)
        iids0 = [i[0] for i in iids0]
        iid_df = pd.DataFrame({"sample": iids0})

        meta_df = pd.merge(
            iid_df, meta_df[["sample", "pop", "super_pop"]], how="left", on="sample")

        assert(len(meta_df) == len(iids))  # Sanity Cehck
        print(f"Loaded {np.shape(meta_df)[0]} individual meta file.")

        return meta_df

    #####################################################
    # Methods to produce artificial individuals

    def create_chunked_roh_individual(self, chunk_length=0.005, roh_list=[], pop_list=["TSI"]):
        """Create a chunked Individual
        chunk_length: Chunk Length [in Morgan].
        roh_list: nx2 table: Column1 & 2: ROH Start and End.
        Return Genotype and list of indices which Individual was copied from"""
        f = self.f
        meta_df = self.meta_df
        iids = self.give_iids(meta_df, pop_list)
        l = len(iids)   # Nr of Fitting Individuals

        print(f'the chromosome to simulate is: {self.ch}')
        if self.ch in ['X', 'x']:
            gts_new, rec = self.create_chunked_gts_X(chunk_length, iids)
        else:
            gts_new, rec = self.create_chunked_gts(chunk_length, iids)
        copy_inds = np.random.randint(l, size=len(roh_list))
        copy_ids = iids[copy_inds]  # Indexes of the Individuals to copy from

        # Random Indiviual to copy from:
        for i in range(len(roh_list)):
            roh_begin, roh_end = roh_list[i]
            id_copy = copy_ids[i]
            gts_new = self.copy_in_roh(gts_new, roh_begin, roh_end, id_copy)

        # Extract names of Copy IIDs
        copy_iids = meta_df["sample"].values[copy_ids]

        return gts_new, copy_iids

    def create_chunked_gts_X(self, chunk_length, iids):
        # for X chromosome simulation only
        f = self.f
        nloci, nind, _ = np.shape(f["calldata/GT"])
        gts_ref = np.array(f["calldata/GT"]).reshape((nloci, 2*nind))

        # load recombination map in morgan
        rec = np.array(f["variants/MAP"]).astype("float")
        ch_min, ch_max = np.min(rec) - 1e-10, np.max(rec)
        nr_chunks = np.ceil((ch_max-ch_min)/chunk_length).astype("int")

        # create all copying indices
        # identify columns that are not the second copy of maleX, which is all missing data
        notMissing = np.intersect1d(np.where(np.min(gts_ref, axis=0) != -1)[0], iids)
        #print(f'not missing: {notMissing}')
        #print(f'number of non-missing haplotype: {len(notMissing)}')
        copy_id = np.random.choice(notMissing, size=2*nr_chunks)
        # Create Length Bin Vector
        len_bins = np.arange(ch_min, ch_max, chunk_length)
        len_bins = np.append(len_bins, ch_max)  # Append the last Value

        if self.output == True:
            print("Setting new Genotypes...")

        # Initialize with invalid Value
        gts_new = -np.ones((nloci, 2), dtype="int")

        for i in range(len(len_bins) - 1):  # Iterate over all Length Bins
            c_min, c_max = len_bins[i], len_bins[i + 1]

            i_min, i_max = np.searchsorted(rec, [c_min, c_max])
            hap1, hap2 = copy_id[2*i], copy_id[2*i+1]
            gts_new[i_min:i_max+1, 1] = gts_ref[i_min:i_max+1, hap1]
            gts_new[i_min:i_max+1, 0] = gts_ref[i_min:i_max+1, hap2]

        if self.output == True:
            print("Finished chunked Genotypes")

        assert(nloci == len(rec))  # Sanity Check
        assert(np.min(gts_new) > -1)  # Sanity Check
        return gts_new, rec


    def create_chunked_gts(self, chunk_length, iids):
        """Create Chunked indiviual from Genotype Matrix f.
        Chunk_Lentgh: Length of the chunking
        iids: The Individual IDs"""
        f = self.f

        k = len(iids)  # Nr of Individuals to copy from
        nr_loci, _, _ = np.shape(f["calldata/GT"])

        # Load Recombination Map in Morgan
        rec = np.array(f["variants/MAP"]).astype("float")
        ch_min, ch_max = np.min(rec) - 1e-10, np.max(rec)

        nr_chunks = np.ceil((ch_max - ch_min) / chunk_length).astype("int")
        # Create all copying indices
        copy_id = np.random.randint(k, size=nr_chunks)
        copy_id = iids[copy_id]  # Extract the IIDs to copy from

        # Create Length Bin Vector
        len_bins = np.arange(ch_min, ch_max, chunk_length)
        len_bins = np.append(len_bins, ch_max)  # Append the last Value
        #print(f'len_bins: {len_bins[-10:]}')

        if self.output == True:
            print("Setting new Genotypes...")

        # Initialize with invalid Value
        gts_new = -np.ones((nr_loci, 2), dtype="int")

        for i in range(len(len_bins) - 1):  # Iterate over all Length Bins
            c_min, c_max = len_bins[i], len_bins[i + 1]

            i_min, i_max = np.searchsorted(rec, [c_min, c_max])
            #print(f'{(i_min, i_max)}:{(c_min, c_max)}')
            ind = copy_id[i]
            gts_new[i_min:i_max + 1,
                    1] = f["calldata/GT"][i_min:i_max + 1, ind, 0]
            gts_new[i_min:i_max + 1,
                    0] = f["calldata/GT"][i_min:i_max + 1, ind, 1]
        # append the last bit if necessary
        if i_max+1 < nr_loci:
            ind = np.random.randint(k)
            gts_new[i_max+1:, 1] = f["calldata/GT"][i_max+1:, ind, 0]
            gts_new[i_max+1:, 0] = f["calldata/GT"][i_max+1:, ind, 1]

        if self.output == True:
            print("Finished chunked Genotypes")

        assert(nr_loci == len(rec))  # Sanity Check
        print(f'missing data in gts_new: {np.where(gts_new<0)}')
        assert(np.min(gts_new) > -1)  # Sanity Check
        return gts_new, rec

    def copy_in_roh(self, gts, roh_begin, roh_end, id_copy):
        """Copy In ROH Block.
        f: HDF where to copy from
        gts: 2xn Matrix of Genotypes to copy ROH block in
        id_copy: Integer Index of which Individual to copy from
        roh_begin, roh_end: Start and End of Block to copy from [in Morgan].
        Return: Genotype Matrix"""
        f = self.f

        # Load Recombination Map in Morgan
        rec = np.array(f["variants/MAP"]).astype("float")
        # Get the Indices where to begin and end ROH
        i_min, i_max = np.searchsorted(rec, [roh_begin, roh_end])

        if self.output == True:
            print(f"Copying in Block: {rec[i_min]:.4f}-{rec[i_max]:.4f} M")

        assert(np.shape(gts)[0] == len(rec))  # Sanity Check

        if np.min(f["calldata/GT"][i_min:i_max+1, id_copy, 1]) == -1:
            gts_copy = f["calldata/GT"][i_min:i_max + 1, id_copy, 0]  # The Stretch to copy in
        else:
            assert(np.min(f["calldata/GT"][i_min:i_max+1, id_copy, 1]) != -1)
            gts_copy = f["calldata/GT"][i_min:i_max + 1, id_copy, 1]  # The Stretch to copy in

        # introduce copying errors
        errWhere = np.where(np.random.rand(i_max + 1 - i_min) <= self.e_rate_ref)[0]
        gts_copy[errWhere] = np.abs(1 - gts_copy[errWhere])
        if len(errWhere)>0:
            print(f'adding {len(errWhere)} errors when copying')

        gts[i_min:i_max + 1, :] = gts_copy[:, None]  # Copy in the Stretch

        assert(np.min(gts) > -1)
        return gts

    def give_iids(self, meta_df="", pop_list=["TSI"]):
        """Return list of Indices in pop_list"""
        if len(meta_df)==0:
            meta_df = self.meta_df

        iids = np.where(meta_df["pop"].isin(pop_list))[0]
        if len(iids) == 0:
            iids = np.where(meta_df["super_pop"].isin(pop_list))[0]
        
        if len(iids) == 0:
            print(f"No individuals in population {pop_list} is found. Please check your input!")
            sys.exit()

        if self.output == True:
            print(f"Found {len(iids)} Individuals in {pop_list}")
        return iids

    def give_popfreq_by_pop(self, conPop):
        """Return allele frequency of the specified population."""
        f = self.f
        if len(conPop) != 0:
            iids = self.give_iids("", conPop) # Just use the existing meta_df, so empty string here
            gts = f["calldata/GT"][:, iids, :] # extract genotypes of relevant individuals
        else:
            gts = np.array(f["calldata/GT"])
        nloci, nind, _ = gts.shape
        gts = gts.reshape(nloci, nind*2)
        # exclude haplotypes with missing genotypes (for male x chromosome, one column is always -1, for example)
        gts = gts.astype(np.float32)
        gts[gts == -1] = np.nan
        return np.nanmean(gts, axis=1)



    def get_gts_pop(self, pop_list, meta_df=""):
        """Find all Individuals from a Population.
        Return their genotypes and iids"""
        f = self.f
        if len(meta_df)==0:
            meta_df = self.meta_df

        iids = np.where(meta_df["pop"].isin(pop_list))[0] # Find individuals

        if self.output == True:
            print(f"Found {len(iids)} Individuals in {pop_list}")

        if len(iids)==0:
            raise RuntimeError(f"No Individuals of {pop_list} Found in Reference")

        gts =  f["calldata/GT"][:, iids, :] # Extract the genotypes
        iids = np.array(f["samples"])[iids]
        return gts, iids


    def produce_allele_counts(self):
        """Creates Allele Counts"""
        raise NotImplementedError()


#########################################
if __name__ == "__main__":
    t = Mosaic_1000G()

    roh_list = [[0.5, 0.55], [0.6, 0.65]]

    gts_roh, copy_iids = t.create_chunked_roh_individual(roh_list=roh_list)

    print(np.shape(gts_roh))
    print(gts_roh[:10])
    print(copy_iids)
    print(len(copy_iids))

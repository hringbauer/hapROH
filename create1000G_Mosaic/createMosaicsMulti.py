"""
Class for preparing 1000 Genome Mosaic Data
Wrapper for createMosaic - to simulate several indiviuals with several ROH
Save the created Data to HDF 5, and the ROH block info to a csv in the same folder.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

from ast import get_docstring
import h5py                            # For Processing HDF5s
import numpy as np
from numpy.core.numeric import _isclose_dispatcher
import pandas as pd
import os                              # To delete File
import sys

from createMosaic import Mosaic_1000G  # Import the


class Mosaic_1000G_Multi(object):
    """Class for Preparing Multiple 1000G Mosaics. And Saving them.
    A wrapper for createMosaic"""

    # Important Parameters:
    ch = 3                              # Which Chromosome to analyze
    pop_list = ["TSI"]                  # Which pop to copy from ["TSI"]
    chunk_length = 0.005                # Chunk Length of Chromosomes 0.005
    roh_lengths = np.ones(5) * 0.05     # ROH Lengths per Individual
    iid = "iid"                         # Prefix of Artificial Individual Names
    n = 3       # Nr of individuals to simulate

    path1000G = "./Data/1000Genomes/HDF5/1240kHDF5/Eur1240chr"
    pop_path = "/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv"
    # Where to save the new HDF5 to by default
    save_path = "./Simulated/1000G_Mosaic/TSI/ch3_5cm/"

    output = True  # whether to Print Output
    m_object = 0  # The Mosaic object

    # parameters for faking readcount data
    cov = 0.0 # genome wide coverage
    con = 0.0 # contamination rate
    err_rate = 1e-3 # sequencing err rate

    def __init__(self, cov=0.0, con=0.0, err_rate=1e-3):
        """Initialize"""
        self.cov = cov
        self.con = con
        self.err_rate = err_rate

    def load_m_object(self):
        """Load the Mosaic Object"""
        # Create Save folder if needed
        print(self.save_path)  # DEBUg

        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        self.m_object = Mosaic_1000G(ch=self.ch, path1000G=self.path1000G,
                                     pop_path=self.pop_path, save_path=self.save_path)

    def save_hdf5(self, gt, ad, ref, alt, pos, rec, samples, path):
        """Create a new HDF5 File with Input Data.
        gt: Genotype data [l,k,2]
        ad: Allele depth [l,k,2]
        ref: Reference Allele [l]
        alt: Alternate Allele [l]
        pos: Position  [l]
        m: Map position [l]
        samples: Sample IDs [k]"""

        l, k, _ = np.shape(gt)  # Nr loci and Nr of Individuals

        if os.path.exists(path):  # Do a Deletion of existing File there
            os.remove(path)

        dt = h5py.special_dtype(vlen=str)  # To have no problem with saving

        with h5py.File(path, 'w') as f0:
            # Create all the Groups
            f_map = f0.create_dataset("variants/MAP", (l,), dtype='f')
            f_ad = f0.create_dataset("calldata/AD", (l, k, 2), dtype='i')
            f_ref = f0.create_dataset("variants/REF", (l,), dtype=dt)
            f_alt = f0.create_dataset("variants/ALT", (l,), dtype=dt)
            f_pos = f0.create_dataset("variants/POS", (l,), dtype='i')
            f_gt = f0.create_dataset("calldata/GT", (l, k, 2), dtype='i')
            f_samples = f0.create_dataset("samples", (k,), dtype=dt)

            # Save the Data
            f_map[:] = rec
            f_ad[:] = ad
            f_ref[:] = ref.astype("S1")
            f_alt[:] = alt.astype("S1")
            f_pos[:] = pos
            f_gt[:] = gt
            f_samples[:] = np.array(samples).astype("S10")

        if self.output == True:
            print(f"Successfully saved {k} individuals to: {path}")

    def create_individuals(self):
        """Creates and saves Genotypes of a number of Individuals.
        roh_lengths: Lengths of the ROH to copy in [LIST]
        Save Individuals as well as ROH positions (in Pandas Dataframe)"""
        m = self.m_object
        chunk_length = self.chunk_length
        pop_list = self.pop_list
        n = self.n
        iid = self.iid
        roh_lengths = self.roh_lengths

        l = len(m.f["variants/MAP"])
        c_min, c_max = np.min(m.f["variants/MAP"]), np.max(m.f["variants/MAP"])
        
        ########################## for debugging purpose ######################
        print(f'c_min: {c_min}, c_max: {c_max}')


        ##########################  end #######################################


        gts = -np.ones((l, n, 2), dtype="int")

        iids = [iid + str(i) for i in range(n)]  # Create IIDs

        # Make the Vectors for Indiviuals
        roh_begins, roh_ends = [], []
        iid_list, copy_iids = [], []

        for i in range(n):  # Iterate over all individuals
            if self.output == True:
                print(f"\nDoing Individual {iids[i]}")

        if self.ch in ['X', 'x']:
            # for male x, we ignore the user input (the number of ROH blocks and the length of each block)
            roh_list = np.array([[c_min, c_max]]) # male X should be ROH everywhere
        else:
            roh_list = self.create_roh_list(roh_lengths, c_min, c_max)
        
        print(f'###### roh length: {roh_lengths} ####')
        print(f'##### roh list: {roh_list} #####')
            
        gts_roh, copy_ids = m.create_chunked_roh_individual(
            chunk_length=chunk_length, roh_list=roh_list, pop_list=pop_list)

        gts[:, i, :] = gts_roh

        # Append ROH information to save
        roh_begins += list(roh_list[:, 0])
        roh_ends += list(roh_list[:, 1])
        iid_list += [iids[i], ] * len(roh_list)
        copy_iids += list(copy_ids)

        # Save the ROH Information:
        if self.output == True:
            print("Finished Copying")
            print("Minimum of Genotypes: ")
            print(np.min(gts))
            print(np.shape(gts))

        # Save the Data:
        self.save_rohlist(roh_begins, roh_ends, iid_list,
                          copy_iids, ch=self.ch)

        self.save_genotypes(m.f, gts, iids)

    def extract_individuals(self):
        """Extracts and saves Genotypes as new hdf5 of a number of Individuals,
        with exactly the same genotypes
        identical state
        """
        m = self.m_object
        pop_list = self.pop_list

        gts, iids = m.get_gts_pop(pop_list)

        self.save_genotypes(m.f, gts, iids)  # Save the genotypes

    def create_roh_list(self, roh_lengths, min_ch, max_ch):
        """Create List of ROH.
        Evenly spaces Chromosome and places ROH randomly in there
        (enforcing gaps at least length of rohs)
        roh_lengths: How long the ROHs are [In Morgan]
        min_ch: Minimum Map Position
        max_ch: Maximum Max Position"""
        k = len(roh_lengths)

        if k == 0:   # If no blocks, return empty list
            return np.array([[], []]).T

        l = (max_ch - min_ch) / k  # Length Batch

        # Factor 2: Make sure that there is some gap
        pos = np.random.random(size=k) * (l - roh_lengths * 2)
        roh_begin = pos + np.arange(k) * l + min_ch
        roh_end = roh_begin + roh_lengths

        roh_list = np.array([roh_begin, roh_end]).T

        return roh_list

    def save_genotypes(self, f, gts, samples, path=""):
        """Save the full genotype Matrix
        f: HDF5 File with Matching entries for Meta Data per Locus
        gts: Genotype Matrix [l,k,2]
        samples: List of sample IIDs [k]
        Path: Where to save to
        cov: genome wide coverage
        con: contamination rate
        err: sequencing error rate"""

        if len(path) == 0:
            path = self.save_path

        path = self.save_path + "data.h5"

        gt = gts
        ad = gts

        ref, alt = f["variants/REF"][:], f["variants/ALT"][:, 0]
        pos = f["variants/POS"]
        rec = f["variants/MAP"]
        print(f'shape of gt: {gt.shape}')

        # add (artificial) readcount coverage data
        cov = self.cov
        con = self.con
        err_rate = self.err_rate

        gt_ref = np.array(self.m_object.f["calldata/GT"])
        nloci_ref, nind_ref, _ = gt_ref.shape
        gt_ref = gt_ref.reshape((nloci_ref, 2*nind_ref))
        assert(nloci_ref == gt.shape[0])
        p = np.mean(gt_ref, axis=1)

        if cov > 0:
            ad = np.zeros_like(gt)
            nloci, nind, _ = ad.shape
            nreads = np.random.poisson(cov, (nloci, nind)) # using Poisson random number to model the number of reads that cover each site
            for i in range(nloci):
                allele_freq = p[i]
                for j in range(nind):
                    gt_ij = gt[i, j]
                    nread_ij = nreads[i, j]
                    # sample reads from 1) the genotype or 2) the reference panel based on contamination rate
                    for read in range(nread_ij):
                        isCon = np.random.rand() <= con
                        if not isCon:
                            # sample reads based on the genotype
                            ret = np.random.choice(gt_ij)
                            if np.random.rand() <= err_rate:
                                ret = abs(1 - ret) # simulate sequencing error
                            ad[i, j, ret] += 1
                        else:
                            # sample reads based on the reference panel
                            ret == 1 if np.random.rand() <= allele_freq else 0
                            if np.random.rand() <= err_rate:
                                ret = abs(1 - ret)
                            ad[i, j, ret] += 1


        # Maybe filter for Loci here
        self.save_hdf5(gt, ad, ref, alt, pos, rec, samples, path)

    def save_rohlist(self, roh_beg, roh_end, iids, copyiids, ch=0, path=""):
        """Save the ROH List to Path"""
        if len(path) == 0:
            path = self.save_path

        if ch == 0:
            ch = self.ch

        # Create the Saving Dataframe
        df_save = pd.DataFrame({"ROH_Begin": roh_beg,
                                "ROH_End": roh_end,
                                "iid": iids,
                                "copyiid": copyiids,
                                "chr": ch})

        path = path + "roh_info.csv"
        df_save.to_csv(path, sep="\t", index=False)

        if self.output == True:
            print(f"Successfully saved to {path}")

    def save_metafile(self, path=""):
        """Save the population File"""
        # if len(path)==0:
        #    path = self.path

        save_df = self.m_object.meta_df
        save_df.to_csv(path, sep="\t", index=False)
        print(f"Successfully saved to {path}")


def multi_run_lengths(base_path="./Simulated/1000G_Mosaic/TSI/", pop_list=["TSI"], n=20,
                      ch=3, chunk_length=0.005, lengths=[1.0, 3.0, 5.0, 10.0], n_blocks=5):
    """Create Multiple ROH runs and saves combined data into base_path hdf5 and roh_info df
    base_path:  Start of SavePaths
    pop_list: The Reference Populations for Mosaic
    n: Number of Individuals to simulate
    chunk_length: Lenths of the Chunks to mosaic
    ch: Chromosome to use
    lengths: The Length of the Blocks to copy in
    n_blocks: The NR of the Blocks to copy in"""

    t = Mosaic_1000G_Multi()  # Create the MultiRUn Object
    t.pop_list = pop_list
    t.n = n
    t.chunk_length = chunk_length
    t.ch = ch  # The Chromosome

    for l in lengths:
        t.roh_lengths = np.ones(n_blocks) * 0.01 * l  # Set the Lengths
        # Where to save the new HDF5 to
        t.save_path = base_path + "ch" + str(t.ch) + "_" + str(int(l)) + "cm/"

        t.load_m_object()
        t.create_individuals()


def copy_population(base_path="./Simulated/1000G_Mosaic/TSI0/",
                    path1000G="./Data/1000Genomes/HDF5/1240kHDF5/Eur1240chr",
                    pop_list=["TSI", ], ch=3):
    """Extract one or more population(s)
    base_path: Where to save extracted data to
    path1000G: Where to take the reference hdf5 from
    pop_list: The copied Populations"""
    t = Mosaic_1000G_Multi()  # Create the MultiRUn Object
    t.pop_list = pop_list
    t.path1000G = path1000G
    t.ch = ch
    t.save_path = base_path + "ch" + str(t.ch) + "/"
    t.load_m_object()  # Load the inner Mosaic Object
    t.extract_individuals()




#############################################
#############################################
# some helper functions to aid simulation
# helper function for the simulation

def create_individual_mosaic(base_path="./Simulated/1000G_Mosaic/TSI/", 
                             path1000G="./Data/1000Genomes/HDF5/1240kHDF5/Eur1240chr",
                    pop_list=["TSI"], n=2, ch=3, chunk_length=0.005, l = 1, n_blocks=5,
                    cov=0.0, con=0.0, err_rate=1e-3):
    """Create Multiple ROH runs and saves combined data into base_path hdf5 and roh_info df
    base_path:  Start of SavePaths
    path1000G: Where to find the 1000 Genome Data
    pop_list: The Reference Populations for Mosaic
    n: Number of Individuals to simulate
    chunk_length: Lenths of the Chunks to mosaic
    ch: Chromosome to use
    l: Length of the ROH blocks
    n_blocks: The NR of the Blocks to copy in
    cov: genome wide coverage
    con: contamination rate
    err_rate: sequencing err rate"""
    
    ########### Pipe the output
    save_path = base_path + "chr" + str(ch) + "_" + str(int(l)) + "cM/"
    
    if not os.path.exists(save_path):
            os.makedirs(save_path)
    
    print(f"Setting save path...: {save_path}")
    sys.stdout = open(save_path + "mosaic_out.txt", 'w')
    
    t = Mosaic_1000G_Multi(cov, con, err_rate)  # Create the MltiRUn Object
    
    ##################################
    ### Set the parameters for the run
    t.pop_list = pop_list
    t.path1000G = path1000G
    t.n = n
    t.chunk_length = chunk_length
    t.ch = ch  # The Chromosome
    t.roh_lengths = np.ones(n_blocks) * 0.01 * l  # Set the Lengths
    t.save_path = save_path
    t.load_m_object()
    t.create_individuals()

def multi_run(fun, prms, processes = 4):
    """Implementation of running in Parallel.
    fun: Function
    prms: The Parameter Files
    processes: How many Processes to use"""
    print(f"Running {len(prms)} jobs in parallel.")
    
    with mp.Pool(processes = processes) as pool:
        results = pool.starmap(fun, prms)










#########################################
if __name__ == "__main__":
    #t = Mosaic_1000G_Multi()
    # t.create_individuals()

    # Do the Multi_Run of Lengths
    # multi_run_lengths()

    ### Test Running FP Individuals without copied blocks:
    #multi_run_lengths(base_path="./Simulated/1000G_Mosaic/TSI4/",
    #                  lengths=[1, ], n_blocks=5, n=2, chunk_length=0.0025)
    
    ### Simulate 100 CHB Mosaic Individuals without copied Blocks
    multi_run_lengths(base_path="./Simulated/1000G_Mosaic/CHB/",
                      lengths=[0], n_blocks=0, n=100, chunk_length=0.0025)

    ### Creating copies of full individuals
    #copy_population(base_path="./Simulated/1000G_Mosaic/CHB0/",
    #                path1000G="./Data/1000Genomes/HDF5/1240kHDF5/NonEur1240chr",
    #                pop_list=["CHB", ], ch=3)

    #copy_population(base_path="./Simulated/1000G_Mosaic/CLM0/",
    #                path1000G="./Data/1000Genomes/HDF5/1240kHDF5/NonEur1240chr",
    #                pop_list=["CLM", ], ch=3)

    #copy_population(base_path="./Simulated/1000G_Mosaic/YRI0/",
    #                path1000G="./Data/1000Genomes/HDF5/1240kHDF5/NonEur1240chr",
    #                pop_list=["YRI", ], ch=3)

"""
Class for loading Eigenstrat Genotype Data from geno file.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""
import numpy as np
import pandas as pd

class EigenstratLoad(object):
    """Class that loads and postprocesses Eigenstrats"""
    base_path = "./Data/ReichLabEigenstrat/Raw/v37.2.1240K"
    nsnp = 0
    nind = 0
    rlen = 0
    output = True
    df_snp = 0  # Dataframe with all SNP

    def __init__(self, base_path="", output=True, sep=r"\s+"):
        """Concstructor:
        base_path: Data Path without .geno/.snp/.ind).
        ch: Which chromosome to load
        sep: What separator to use when loading the File"""
        self.output = output
        if len(base_path) > 0:
            self.base_path = base_path
        geno_file = open(self.base_path + ".geno", "rb")
        header = geno_file.read(20)  # Ignoring hashes
        self.nind, self.nsnp = [int(x) for x in header.split()[1:3]]
        # assuming sizeof(char)=1 here
        self.rlen = max(48, int(np.ceil(self.nind * 2 / 8)))

        self.df_snp = self.load_snp_df(sep=sep)   # Load the SNP DataFrame
        self.df_ind = self.load_ind_df(sep=sep)   # Load the Individual DataFrame
        assert(len(self.df_snp) == self.nsnp)  # Sanity Check
        assert(len(self.df_ind) == self.nind)  # Sanity Check II

        if self.output == True:
            print(f"3 Eigenstrat Files with {self.nind} Individuals and {self.nsnp} SNPs")

    def load_snp_df(self, sep=r"\s+"):
        """Load the SNP dataframe.
        Uses self.base_path
        sep: What separator to use when loading the File"""
        path_snp = self.base_path + ".snp"
        df_snp = pd.read_csv(path_snp, header=None,
                             sep=sep, engine="python")
        df_snp.columns = ["SNP", "chr", "map",
                          "pos", "ref", "alt"]  # Set the Columns
        return df_snp

    def load_ind_df(self, sep=r"\s+"):
        """Load the Individual dataframe.
        Uses self.base_path
        sep: What separator to use when loading the File"""
        path_ind = self.base_path + ".ind"
        df_ind = pd.read_csv(path_ind, header=None,
                             sep=r"\s+", engine="python")
        df_ind.columns = ["iid", "sex", "cls"]  # Set the Columns
        return df_ind
    
    def get_geno_all(self, missing_val=3):
        """Load all genotypes from Eigenstrat File.
        Use self.nind for number of individuals.
        Return genotype matrix, with missing values set to missing_val"""
        geno = self.give_bit_file()  # Load the whole bit file
        gt = np.unpackbits(geno, axis=1)[:,:2*self.nind]
        gt = 2 * gt[:, 0::2] + gt[:, 1::2]
        gt[gt == 3] = missing_val  # set missing values
        return gt

    def get_geno_i(self, i, missing_val=3):
        """Load Individual i"""
        batch, eff = self.get_enc_index(i)
        geno = self.give_bit_file()  # Load the whole bit file

        geno_sub = geno[:, [batch]]  # Byte value of batch
        geno_sub = np.unpackbits(geno_sub, axis=1)[:, 2 * eff:2 * eff + 2]
        geno_sub = 2 * geno_sub[:, 0] + geno_sub[:, 1]
        geno_sub[geno_sub == 3] = missing_val  # set missing values
        return geno_sub

    def get_geno_indvidiaul_iid(self, iid):
        """Return Genotypes of Individual iid"""
        i = self.get_index_iid(iid)
        g = self.get_geno_i(i)
        return g

    def give_bit_file(self):
        base_path = self.base_path
        geno = np.fromfile(self.base_path + ".geno",
                           dtype='uint8')[self.rlen:]  # without header
        geno.shape = (self.nsnp, self.rlen)
        return geno

    def get_enc_index(self, i):
        """Get the Index in the Encoding and the modulo 4 value
        (position in batch)"""
        rlen_sub = int(np.floor(i * 2 / 8))  # Effectively dividing by four
        mod_i = i % 4  # Calculate the rest
        return rlen_sub, mod_i

    def give_positions(self, ch):
        """Return Array of Positions and Indices of
         all SNPs on Chromosome ch"""
        df_snp = self.df_snp
        ch_loci = (df_snp["chr"] == ch)
        idcs = np.where(ch_loci)[0]
        assert(len(idcs) > 0)
        pos = df_snp.loc[ch_loci, "pos"]
        return pos, idcs

    def give_ref_alt(self, ch):
        """Return Arrays of Ref/Alt of all SNPs on Chromosome ch"""
        df_snp = self.df_snp
        df_t = df_snp.loc[df_snp["chr"] == ch, ["ref", "alt"]]  # Subset to ch
        ref, alt = df_t["ref"].values, df_t["alt"].values
        return ref, alt

    def get_index_iid(self, iid):
        """Get Index of Individual iid"""
        # Detect the Individual
        found = np.where(self.df_ind["iid"] == iid)[0]
        if len(found)==0:
            raise RuntimeError(f"Individual {iid} not found!")
        else: 
            i = found[0]
        return i

    def extract_snps(self, id, markers, conversion=True, dtype=np.int8):
        """Extract SNPs for Integer Index i on marker list
        markers. If conversion: Convert to VCF type encoding
        Load all SNPs and then subset"""
        geno = self.get_geno_i(id)
        geno = geno[markers] # Subset to markers

        if conversion == True:
            geno_new = -np.ones((2,len(geno)), dtype=dtype)
            geno_new[:, geno==0]=1    # 2 Derived Alleles
            geno_new[:, geno==2]=0    # 2 Ancestral Alleles
            ### Heterozgyotes
            geno_new[0, geno==1]=1
            geno_new[1, geno==1]=0
            geno = geno_new

        return geno

#########################################################
#########################################################
#### Subclass for Non-Binary Eigenstrats

class EigenstratLoadUnpacked(EigenstratLoad):
    """Class that loads and postprocesses Eigenstrats.
    Same as Superclass, but overwrites methods to load
    non-binary encoded Genotype Data"""

    def __init__(self, base_path="", output=True, sep=r"\s+"):
        """Overwrite Concstructor:
        base_path: Data path without .geno/.snp/.ind.
        ch: Which chromosome to load
        sep: What separator to use when loading the File"""
        self.output = output
        if len(base_path) > 0:
            self.base_path = base_path
        ### Get Size of Data Matrix and sanity check
        with open(self.base_path + ".geno",'r') as f:
            t = f.read()
            l = t.splitlines()
            self.nind = len(l[0])
            self.nsnp = len(l)

        self.df_snp = self.load_snp_df(sep=sep)   # Load the SNP DataFrame
        self.df_ind = self.load_ind_df(sep=sep)   # Load the Individual DataFrame
        assert(len(self.df_snp) == self.nsnp)  # Sanity Check
        assert(len(self.df_ind) == self.nind)  # Sanity Check II

        if self.output == True:
            print(f"3 Eigenstrat Files with {self.nind} Individuals and {self.nsnp} SNPs")

    def get_geno_i(self, i, missing_val=3):
        """Load Genotype for Individual (Row) i,
        assuming it's encoded in unpacked Format."""
        geno=np.genfromtxt(self.base_path + ".geno", dtype='i1', delimiter=1, usecols=i)
        geno[geno == 9] = missing_val
        return geno

#########################################################
#########################################################
def is_binary_file(path, extension=".geno"):
    """Test whether a file at path + extension is binary.
    Return boolean if the case """
    binary=False
    try:
        with open(path + extension, "r") as f:
            t = f.readline()
    except UnicodeDecodeError:
        binary=True
    return binary

#########################################################
#########################################################
#### Factory Method

def load_eigenstrat(base_path, output=True, sep=r"\s+", packed=-1):
    """Factory Method to Load Eigenstrat object
    sep: What separator to use when loading the File. 
    Default is space seperated (by arbitrary number of spaces)
    Packed: Whether Genotype Data is encoded in binary Format"""
    
    ### Determine automatically
    if packed==-1:
        packed = is_binary_file(base_path, extension=".geno")
        if output:
            print(f"Eigenstrat packed: {packed}")
    
    ### Load
    if packed:
    	es = EigenstratLoad(base_path, output=output, sep=sep)
    else:
    	es = EigenstratLoadUnpacked(base_path, output=output, sep=sep) 
    return es


#########################################################
#########################################################
# Some Testing
if __name__ == "__main__":
    es = load_eigenstrat(base_path="./Data/ReichLabEigenstrat/Olalde2019/Olalde_et_al_genotypes")
    #i = es.get_index_iid(iid="I4055")
    pos, idcs = es.give_positions(ch=3)
    print(idcs)
    print(len(pos))
    print(len(idcs))

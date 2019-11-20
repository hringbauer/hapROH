import pandas as pd
import numpy as np

class Expected_Roh():
    """Class that calculates expected ROH"""
    # Length of all Chromosomes
    chr_lgts = [2.8426, 2.688187, 2.232549, 2.14201, 2.040477, 1.917145, 1.871491, 1.680018, 
                1.661367, 1.8090949, 1.5821669, 1.745901, 1.2551429, 1.1859521, 1.413411, 
                1.340264, 1.2849959, 1.175495, 1.0772971, 1.082123, 0.636394, 0.724438]

#################################    
#################################
### For relatives m Meiosis apart

    def roh_pdf_rel(self, x, chr_l, m):
        """Gives back the pdfs for Blocks of Length l [Morgan]
        on a Chromosome of Length [Morgan].
        m: Nr of Meiosis.
        Return PDF (per Morgan)"""
        pdf0 = (chr_l - x) * m**2 * np.exp(-x * m) # In Chromosome
        pdf1 = 2 * np.exp(- x * m) * m    # Boundary
        return (pdf0 + pdf1) * (x < chr_l)  # If x < chr_l return 0

    def roh_pdf_allchr_rel(self, x, chr_lgts, m):
        """Calculate the PDF of ROH blocks of length x [Morgan]
        for m Recombination events summed over Chromosome Lengths
        x: Can be Array
        chr_lgts: Array of all chromosome lengths [in Morgan]
        Return PDF (per Morgan)"""
        pdfs = [self.roh_pdf_rel(x, chr_l, m) for chr_l in chr_lgts]
        pdf_full = np.sum(pdfs, axis=0)
        return pdf_full

    def coal_prob_rel(self, m, comm_anc=1):
        """Calculate Coalescence Probability.
        m: Nr of Meiosis
        comm_anc: How many common ancestors"""
        c_prob = comm_anc * (1 / 2) ** m
        return c_prob

    def exp_roh_ind_rel(self, x, m, comm_anc=1, chr_lgts=[]):
        """Calculates the Expected density of ROH Blocks per Morgan at x M
        for full Individual
        x: Array of Block Lenths
        m: Nr of Meisois
        comm_anc: Nr of Ancestry Loops"""  
        if len(chr_lgts)==0:
            chr_lgts = self.chr_lgts ## Default value
        pdf_full = self.roh_pdf_allchr_rel(x, chr_lgts, m)
        c_prob = self.coal_prob_rel(m, comm_anc)
        exp_blocks = pdf_full * c_prob
        return exp_blocks

    def exp_roh_len_in_bin_rel(self, l=[0,3], m=6, comm_anc=4, bins=1000):
        """Calculate Epected Block Length for Individual m Meisosis apart and
        with comm_anc haplotypes in common within interval l [in Morgan].
        bins: Number of linearly spaced bins to do.
        Return Sum [in Morgan]"""
        x_arr = np.linspace(l[0], l[1], bins)
        bw = x_arr[1] - x_arr[0]
        y = self.exp_roh_ind_rel(x_arr + bw/2, m=m, comm_anc=comm_anc) * (x_arr+bw/2) * bw
        return np.sum(y)
    
#################################    
#################################
### For Pop Size N (=2Ne)

    def roh_pdf_N(self, x, chr_l, N):
        """Gives back the pdfs for Blocks of Length l [Morgan]
        on a Chromosome of Length [Morgan].
        m: Nr of Meiosis.
        Return PDF (per Morgan)"""
        pdf0 = 8 * (chr_l - x)/N * 1/(2*x + 1/N)**3  # In Chromosome
        pdf1 = 4/N * (1 / (2*x + 1/N)**2)     # Boundary
        return (pdf0 + pdf1) * (x < chr_l) * (x>0)  # If x < chr_l return 0, no such blocks

    def roh_pdf_allchr_N(self, x, chr_lgts=[], N=1000):
        """Calculate the PDF of ROH blocks of length x [Morgan]
        for m Recombination events summed over Chromosome Lengths
        x: Can be Array
        chr_lgts: Array of all chromosome lengths [in Morgan]
        Return PDF (per Morgan)"""
        if len(chr_lgts)==0:
            chr_lgts = self.chr_lgts ## Default value
        pdfs = [self.roh_pdf_N(x, chr_l, N) for chr_l in chr_lgts]
        pdf_full = np.sum(pdfs, axis=0)
        return pdf_full

    def exp_roh_len_in_bin_N(self, l=[0,3], bins=1000, N=500, chr_lgts=[]):
        """Calculate Epected Block Length for Individual m Meisosis apart and
        with comm_anc haplotypes in common within interval l [in Morgan].
        bins: Number of linearly spaced bins to do.
        N: Population size of haploids
        Return Sum [in Morgan]"""
        if len(chr_lgts)==0:
            chr_lgts = self.chr_lgts ## Default value
        x_arr = np.linspace(l[0], l[1], bins)
        bw = x_arr[1] - x_arr[0]
        y = self.roh_pdf_allchr_N(x_arr + bw/2, chr_lgts=chr_lgts, N=N) * (x_arr+bw/2) * bw
        return np.sum(y)
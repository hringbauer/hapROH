"""
Contains various Functions to operate with h5 files:
Loading h5, as well as converting to VCFs
@ Author: Harald Ringbauer, 2019, All rights reserved
"""

import numpy as np
import pandas as pd
import h5py
from scipy.stats import binom  # Binomial Likelihood

def load_h5(path, output=True):
    """Load HDF5 from path and return hdf5 object"""
    f = h5py.File(path, "r") # Load for Sanity Check. See below!

    if output == True:
        print("Loaded HDF5")
        print("Loaded %i variants" % np.shape(f["calldata/GT"])[0])
        print("Loaded %i individuals" % np.shape(f["calldata/GT"])[1])
        print(list(f["calldata"].keys()))
        print(list(f["variants"].keys()))
        #self.f["samples"] # Samples Vector

    ### Sanity Check whether both Genotypes are there and nothing else
    assert(np.min(f["calldata/GT"]) >= -1)
    assert(np.max(f["calldata/GT"]) == 1)
    
    return f


def to_vcf(chrom, pos, ref, alt, gt, iids, vcf_path, header=[], pl=[]):
    """Saves VCF. If Genotype Likelihoods given (pl), save them too."""
    ### Hard-Coded Default Header
    if len(header)==0:
        header = """##fileformat=VCFv4.3\n##FILTER=<ID=PASS,Description="All filters passed">\n##fileDate=20191010\n##source=1000GenomesPhase3Pipeline\n##reference=ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz\n##contig=<ID=3,assembly=b37,length=198022430>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n"""     
        
    #last_line_h =  "\n#CHROM POS ID REF ALT QUAL FILTER INFO"
    dct = {'#CHROM':chrom, 'POS':pos, 'REF':ref, 'ALT':alt}
    df = pd.DataFrame(dct)
    df['ID'] = ""
    df['QUAL'] = 40
    df['FILTER'] = "PASS"
    df['INFO'] = ""
    
    frmt = "GT"
    if len(pl)>0:
        frmt = "GT:PL"
        header=header + """##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">\n"""
    df["FORMAT"] = frmt

    df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', "FORMAT"]] 

    ### Add the Genotype Data
    add_gt_data(df, gt, pl=pl, iids=iids)
    
    ### Write the Header
    with open(vcf_path, 'w') as vcf:
        vcf.write(header)
        #vcf.write(last_line_h)

    #### Write the tab seperated data
    df.to_csv(vcf_path, sep="\t", mode='a', index=False)  # Append
    print(f"Successfully saved VCF to {vcf_path}")
    
def add_gt_data(df, gt, pl=[], iids=[], m_sym="."):
    """Add Genotype and Allele Depth Fields [l,n,2] for iids to pandas dataframe df.
    Return modified Data Frame".
    If pl (Genotype Likelihoods) given, add them too."""
    assert(np.shape(gt)[1]==len(iids)) # Sanity Check
    
    ### Replace missing Data with dot again
    missing = gt<0  # Missing Data
    gt = gt.astype("str")
    gt[missing] = m_sym
    
    gt_vcf = np.core.defchararray.add(gt[:,:,0], "/")
    gt_vcf = np.core.defchararray.add(gt_vcf, gt[:,:,1])
    
    if len(pl)>0:
        pl = pl.astype("str")
        gt_vcf = np.core.defchararray.add(gt_vcf, ":")
        gt_vcf = np.core.defchararray.add(gt_vcf, pl[:,:,0])
        gt_vcf = np.core.defchararray.add(gt_vcf, ",")
        gt_vcf = np.core.defchararray.add(gt_vcf, pl[:,:,1])
        gt_vcf = np.core.defchararray.add(gt_vcf, ",")
        gt_vcf = np.core.defchararray.add(gt_vcf, pl[:,:,2])
        
    for i, iid in enumerate(iids):
        #data = map('/'.join, zip(gt[:,i,0], gt[:,i,1]))
        df[iid] = gt_vcf[:,i]
        #if len(ad)>0:   # Add Allele Depth Data if needed
        #    print("Implement this") 
    return df
    

def hdf5_to_vcf(path_h5, path_vcf, iids=[], markers=[], chrom=0, pl_field=False):
    """Load HDF5 from path_h5, extract iids and
    (if given) markers by position and save vcf to path_vcf.
    pl: If True, also save Genotype Likelihoods!
    chrom: Value for chromosome (otherwise load from h5)
    iids: Which Individuals to match and save. If none given: Save all!"""
    
    f = load_h5(path=path_h5)
    
    if len(iids)==0:
        iids = f["samples"][:]
        
    if chrom==0:
        chrom = f["variants/CHROM"][:]
        
    pos = f["variants/POS"][:]
    ref = f["variants/REF"][:] 
    alt = f["variants/ALT"][:] 
    
    idx = np.isin(f["samples"], iids)
    gt = f["calldata/GT"][:,idx,:]
    
    ### Get Genotype Likelihoods from AD field
    pl=[]  # Default
    if pl_field:  
        print("Calculating Genotype Likelihoods...")
        ad = f["calldata/AD"][:]
        gl = ad_to_gentoypeL(ad)  # Convert Allele Depths to Genotype Likelihood
        pl = gl_to_pl(gl)         # Convert Genotype Likelihood to PHRED scale 
        
    print("Saving to VCF...")
    to_vcf(chrom, pos, ref, alt, gt, iids, path_vcf, pl=pl)


##################################################
### Functions for Genotype Likelihoods from AD

def ad_to_gentoypeL(ad, error=0.001):
    """Convert Allele Depth Fields to Genotype Likelihoods.
    ad: [l,n,2] contains allele contains readcounts (integers)
    error: Flip Error for Read
    return: Genotype Probabilities (Pr(G|RC)) [l,n,3] for 00/01/11"""
    rc_tot = np.sum(ad, axis=2)
    rc_der = ad[:,:,1]

    p_read = np.array([error, 0.5, 1-error])  # Probabilities for the 3 Genotypes
    prob_binom = binom.pmf(rc_der[:, :, None], rc_tot[:, :, None], p_read[None, None, :])
    return prob_binom

def gl_to_pl(gl):
    """Convert Genotype Probabilities to normalized PHRED scores
    gl: [l,n,3] Probabilities Pr(G|RC) (not logscale)
    return: [l,n,3] vector"""
    gl = -10 * np.log10(gl)  # Convert to PHRED scale
    assert(np.min(gl)>=0)
    pl = gl - np.min(gl, axis=2)[:,:,None] # Normalize
    pl = np.round(pl).astype("int16")  # Round to Integers
    assert(np.min(pl)>=0) # Sanity Check
    pl = np.clip(pl, a_min=0, a_max=99)  # Clip to 99    
    return pl



##################################################################
### Example Case, test whether ad_to_genotypeL does what it should
if __name__ == '__main__':
	ad = np.array([[[1,0], [1,1]],  [[3,3], [3,0]], [[10,10], [0,10]]])

	gl = ad_to_gentoypeL(ad)
	pl = gl_to_pl(gl)

	i,j=2,1
	print(ad[i, j, :])
	print(gl[i, j, :])
	print(pl[i, j, :])








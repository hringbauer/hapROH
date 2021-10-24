"""
Create MAF filtered HDF5 Reference Files from 1000 Genome Autosomes
Called with array jobs from sbatch
@ Author: Harald Ringbauer, 2021
"""

import allel
import h5py  # Python Package to do the HDF5.
import numpy as np
import pandas as pd
import socket
import os as os
import sys as sys


#### 1) Set the Path to default HAPSBURG
path = "/project2/jnovembre/hringbauer/HAPSBURG/"  # The Path on Midway Cluster
os.chdir(path)
print(os.getcwd())
sys.path.insert(0,"./package/")  # hack to get local package first in path
from hapsburg.PackagesSupport.h5_python.h5_functions import merge_in_ld_map, save_data_h5

#########################################################
### Helper Functions

def filter_vcf_to_biallelic(in_path, out_path):
    """Filters .vcf File to biallelic SNPs"""
    command = f"bcftools view -Oz -o {out_path} -m2 -M2 -v snps {in_path}"
    os.system(command)
    #!bcftools view -Oz -o $out_path -m2 -M2 -v snps $in_path 
    print("Finished BCF tools filtering to biallelic variants.")

def vcf_to_hdf5(in_path, out_path, path_vcf100g=""):
    """Transform Full VCF to full HDF5"""
    allel.vcf_to_hdf5(input=in_path, output=out_path, compression="gzip") 
    print("Finished Transforming VCF to HDF5.")
    # Takes 10 Minutes ####chunk_length=1000000, chunk_width=1, garbage performance
    
def download_1kg(path_source="", path_target=""):
    """cluster: Whether program is run on cluster"""
    tbi_s = path_source + ".tbi"
    tbi_t = path_target + ".tbi"
    
    #!scp $path_source $path_target # Only Download the .vcf (not the .tbi)
    #!scp $tbi_s $tbi_t             # Download the tbi
    
    command = f"scp {path_source} {path_target}"
    command1 = f"scp {tbi_s} {tbi_t}" 
    os.system(command)
    os.system(command1)
    
    print(f"Transfer complete. To {path_target}")
    
def downsample_af(path_h5, path_save, maf=0.05, batch=1000000):
    """Loads hdf5 at path_h5, and filters to loci >maf.
    Save h5 at path_save. Assumes everything is in standard format"""

    f = h5py.File(path_h5, "r") # Load for Sanity Check. See below

    ### Calculate the AFs
    k = np.shape(f["calldata/GT"])[0]
    
    idcs = [] # the list of list of idx in each batch
    ### Work in batches
    for i in range(int(k/batch)+1):
        print(f"Doing batch {i}...")
        gt = f["calldata/GT"][i*batch:(i+1)*batch,:,:]  
        assert(np.min(gt)==0) # Sanity Check for old data
        assert(np.max(gt)==1) # Sanity Check for old data
        n = np.shape(gt)[1]*2 # number of haplotypes
        gt = np.sum(gt==0, axis=1)
        gt = np.sum(gt, axis=1)
        p_der = 1 - gt/n
        idx = (p_der > maf) & (p_der < 1 - maf) # Filter for MAF
        idcs.append(idx)
        
    idx = np.concatenate(idcs)
    print(f"Downsampling to {np.sum(idx)}/{len(idx)} Markers with MAF >{maf}")

    gt = f["calldata/GT"][:,:,:][idx] # Extract the Ind
    assert(np.min(gt)==0) # Sanity Check for new data
    assert(np.max(gt)==1) # Sanity Check for new data
    
    ref=f["variants/REF"][:][idx]
    if len(np.shape(f["variants/ALT"]))>1:
        alt=f["variants/ALT"][:][idx,0] # only save the first alt allele
    else:
        alt=f["variants/ALT"][:][idx]
        
    pos=f["variants/POS"][:][idx]
    rec=f["variants/MAP"][:][idx]
    samples=f["samples"][:]
    f.close()
    
    print("Saving new HDF...")
    save_data_h5(gt=gt, ad=[],
                 ref=ref, alt=alt,
                 pos=pos, rec=rec,
                 samples=samples,
                 path=path_save,
                 compression='gzip',
                 ad_group=False, gt_type='int8')
    
#########################################################
### Master Function

def full_prep_h5_full(ch = 4, maf=0.05):
    """Function to run the full preparation of HDF from 1000G VCF. 
    Comment out steps if not needed."""
    ### The Paths
    path_vcf_source = f"/project2/jnovembre/data/external_public/1kg_phase3/haps/ALL.chr{ch}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    path_vcf_here = f"./Data/1000Genomes/AutosomeVCF/Full/chr{ch}.vcf.gz"
    path_vcf_filter = f"./Data/1000Genomes/AutosomeVCF/Full/chr{ch}.biallelic.vcf.gz"
    path_h5 = f"./Data/1000Genomes/HDF5/FULLHDF5/chr{ch}.hdf5"
    path_h5_maf = f"./Data/1000Genomes/HDF5/FULLHDF5/maf5_chr{ch}.hdf5"
    path_snp1240k = "./Data/1000Genomes/Markers/MinMyc.snp"
    
    ### The functions of the pipeline
    download_1kg(path_vcf_source, path_vcf_here)  ## Takes about 20 seconds
    filter_vcf_to_biallelic(path_vcf_here, path_vcf_filter) ### Takes about XX seconds
    vcf_to_hdf5(path_vcf_filter, path_h5)            ## Takes about 20 minutes
    merge_in_ld_map(path_h5=path_h5, path_snp1240k=path_snp1240k, chs=[ch,], write_mode="a")  ## Takes about XX seconds, mode a for new LD Map
    downsample_af(path_h5, path_h5_maf, maf=maf)
    
    print(f"Full Run Finished Successfully.")


#########################################################
#########################################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Script needs argument (indiviual i)")
    ch = int(sys.argv[1]) # The Parameter passed to the Python Script from outside
    
    full_prep_h5_full(ch = ch, maf=0.05)
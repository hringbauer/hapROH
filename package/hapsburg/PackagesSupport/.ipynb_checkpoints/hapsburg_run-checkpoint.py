"""
Function to run HAPSBURG on a full individual and full reference Dataset,
with all relevant keywords.
Function for running on single chromsome, with all relevant keywords.
@Author: Harald Ringbauer, 2019
"""

import numpy as np
import multiprocessing as mp
import pandas as pd

from hapsburg.hmm_inference import HMM_Analyze   # The HMM core object
from hapsburg.PackagesSupport.parallel_runs.helper_functions import prepare_path, multi_run, combine_individual_data


def hapsb_chrom(iid, ch=3, save=True, save_fp=False, n_ref=2504, exclude_pops=[],
                e_model="EigenstratPacked", p_model="MosaicHDF5", readcounts=True, random_allele=True,
                post_model="Standard", path_targets = "./Data/SA_1240kHDF5/IPK12.h5",
                h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/all1240/chr", 
                meta_path_ref = "./Data/1000Genomes/Individuals/meta_df_all.csv",
                folder_out="./Empirical/Eigenstrat/Reichall/test/", prefix_out="",
                roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.0,
                max_gap=0, cutoff_post = 0.999, roh_min_l = 0.01, logfile=True):
    """Run Hapsburg analysis for one chromosome on eigenstrat data
    Wrapper for HMM Class.
    iid: IID of the Target Individual, as found in Eigenstrat [str]
    ch: Chromosome to run [float]
    path_targets: Path of the target files [str]
    h5_path1000g: Path of the reference genotypes [str]
    meta_path_ref: Path of the meta file for the references [str]
    folder_out: Path of the basis folder for output [str]
    prefix_out: Path to insert in output string, e.g. test/ [str]
    e_model: Emission model to use [str]
    p_model: Preprocessing model tu use [str]
    post_model: Model to post-process the data [str]
    processes: How many Processes to use [int]
    delete: Whether to delete raw posterior per locus [bool]
    output: Whether to print extensive output [bool]
    save: Whether to save the inferred ROH [bool]
    save_fp: Whether to save the full posterior matrix [bool]
    n_ref: Number of (diploid) reference Individuals to use [int]
    exclude_pops: Which populations to exclude from reference [list of str]
    readcounts: Whether to load readcount data [bool]
    random_allele: Whether to pick a random of the two target alleles per locus [bool]
    roh_in: Parater to jump into ROH state (per Morgan) [float]
    roh_out: Parameter to jump out of ROH state (per Morgan) [float]
    roh_jump: Parameter to jump (per Morgan) [float]
    e_rate: Error rate target [float]
    e_rate_ref: Error rate refernce [float]
    cutoff_post: Posterior cutoff [float]
    max_gap: Maximum gap to merge (Morgan) [float]
    roh_min_l: Minimum length ROH (Morgan) [float]
    logfile: Whether to use logfile [bool]
    combine: Wether to combine output of all chromosomes [bool]
    file_result: Appendix to individual results [string]"""
    
    ### Create Folder if needed, and pipe output if wanted
    path_out = prepare_path(folder_out, iid, ch, prefix_out, logfile=logfile)
    
    hmm = HMM_Analyze(cython=2, p_model=p_model, e_model=e_model,
                      manual_load=True, save=save, save_fp=save_fp)

    ### Load and prepare the pre-processing Model
    hmm.load_preprocessing_model()              # Load the preprocessing Model
    hmm.p_obj.set_params(readcounts = readcounts, random_allele=random_allele,
                         folder_out=folder_out, prefix_out_data=prefix_out, 
                         excluded=exclude_pops)
    
    ### Set the paths to ref & target
    hmm.p_obj.set_params(h5_path1000g = h5_path1000g, path_targets = path_targets, 
                         meta_path_ref = meta_path_ref, n_ref=n_ref)
    hmm.load_data(iid=iid, ch=ch)  # Load the actual Data
    hmm.load_secondary_objects()

    ### Set the Parameters
    hmm.e_obj.set_params(e_rate = e_rate, e_rate_ref = e_rate_ref)
    hmm.t_obj.set_params(roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump)
    hmm.post_obj.set_params(max_gap=max_gap, cutoff_post=cutoff_post, roh_min_l = roh_min_l)
    
    ### hmm.calc_viterbi_path(save=save)           # Calculate the Viterbi Path.
    hmm.calc_posterior(save=save)              # Calculate the Posterior.
    hmm.post_processing(save=save)             # Do the Post-Processing.
         

#########################################################
### Run Hapsburg for one Individual (wrap for Chr.)

def hapsb_ind(iid, chs=range(1,23), 
              path_targets = "./Data/ReichLabEigenstrat/Raw/v37.2.1240K",
              h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/all1240int8/chr", 
              meta_path_ref = "./Data/1000Genomes/Individuals/meta_df_all.csv",
              folder_out="./Empirical/Eigenstrat/Reichall/test/", prefix_out="",
              e_model="haploid", p_model="EigenstratPacked", post_model="Standard",
              processes=1, delete=False, output=True, save=True, save_fp=False, 
              n_ref=2504, exclude_pops=[], readcounts=True, random_allele=True,
              roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.00, 
              cutoff_post = 0.999, max_gap=0, roh_min_l = 0.01, logfile=True, combine=True, 
              file_result="_roh_full.csv"):
    """Analyze a full single individual in a parallelized fasion. Run all Chromosome analyses in parallel
    Wrapper for hapsb_chrom
    iid: IID of the Target Individual, as found in Eigenstrat [str]
    path_targets: Path of the target files [str]
    h5_path1000g: Path of the reference genotypes [str]
    meta_path_ref: Path of the meta file for the references [str]
    folder_out: Path of the basis folder for output [str]
    prefix_out: Path to insert in output string, e.g. test/ [str]
    e_model: Emission model to use [str]: haploid/diploid_gt/readcount
    p_model: Preprocessing model tu use [str]: EigenstratPacked/EigenstratUnpacked/MosaicHDF5
    post_model: Model to post-process the data [str]: Standard/MMR (experimental)
    processes: How many Processes to use [int]
    delete: Whether to delete raw posterior per locus [bool]
    output: Whether to print extensive output [bool]
    save: Whether to save the inferred ROH [bool]
    save_fp: Whether to save the full posterior matrix [bool]
    n_ref: Number of (diploid) reference Individuals to use [int]
    exclude_pops: Which populations to exclude from reference [list of str]
    readcounts: Whether to load readcount data [bool]
    random_allele: Whether to pick a random of the two target alleles per locus [bool]
    roh_in: Parater to jump into ROH state (per Morgan) [float]
    roh_out: Parameter to jump out of ROH state (per Morgan) [float]
    roh_jump: Parameter to jump (per Morgan) [float]
    e_rate: Error rate target [float]
    e_rate_ref: Error rate refernce [float]
    cutoff_post: Posterior cutoff [float]
    max_gap: Maximum gap to merge (Morgan) [float]
    roh_min_l: Minimum length ROH (Morgan) [float]
    logfile: Whether to use logfile [bool]
    combine: Wether to combine output of all chromosomes [bool]
    file_result: Appendix to individual results [string]
    
    default is with default Parameters finetuned from 1240k data,
    applied to a 1240k Eigenstrat pseudo-haploid dataset."""
                            
    if output:
        print(f"Doing Individual {iid}...")
    
    ### Prepare the Parameters for that Indivdiual
    prms = [[iid, ch, save, save_fp, n_ref, exclude_pops, e_model, p_model, readcounts, random_allele,
            post_model, path_targets, h5_path1000g, meta_path_ref, folder_out, prefix_out,
            roh_in, roh_out, roh_jump, e_rate, e_rate_ref, max_gap, cutoff_post, roh_min_l, logfile] for ch in chs]
    assert(len(prms[0])==25)   # Sanity Check
                            
    ### Run the analysis in parallel
    multi_run(hapsb_chrom, prms, processes = processes)
                            
    ### Merge results for that Individual
    if combine:
        if output:
            print(f"Combining Information for {len(chs)} Chromosomes...")
        combine_individual_data(folder_out, iid=iid, delete=delete, chs=chs, 
                                prefix_out=prefix_out, file_result=file_result)
    if output:
        print(f"Run finished successfully!")
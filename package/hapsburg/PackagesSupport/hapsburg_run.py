"""
Function to run HAPSBURG on a full individual and full reference Dataset,
with all relevant keywords.
Function for running on single chromsome, with all relevant keywords.
@Author: Harald Ringbauer, 2019
"""

import math
import numpy as np
import multiprocessing as mp
import pandas as pd
import sys
from scipy.optimize import minimize
import numdifftools as ndt


from hapsburg.hmm_inference import HMM_Analyze   # The HMM core object
from hapsburg.PackagesSupport.parallel_runs.helper_functions import prepare_path, multi_run, combine_individual_data, move_X_to_parent_folder

def hapsb_chunk_negloglik(iid, ch, start, end, path_targets, h5_path1000g, meta_path_ref,
                folder_out, c, conPop=["CEU"], roh_in=1, roh_out=0, roh_jump=300, e_rate=0.01, e_rate_ref=1e-3,
                save=False, save_fp=False, n_ref=2504, diploid_ref=True, 
                exclude_pops=[], e_model="readcount_contam", p_model="SardHDF5", 
                readcounts=True, random_allele=False, prefix_out="", logfile=False):
    parameters = locals() # Gets dictionary of all local variables at this point
    
    ### Create Folder if needed, and pipe output if wanted
    _ = prepare_path(folder_out, iid, ch, prefix_out, logfile=logfile) # Set the logfile
    hmm = HMM_Analyze(cython=3, p_model=p_model, e_model=e_model, output=False,
                      manual_load=True, save=save, save_fp=save_fp, start=start, end=end)

    ### Load and prepare the pre-processing Model
    hmm.load_preprocessing_model(conPop)              # Load the preprocessing Model
    hmm.p_obj.set_params(readcounts = readcounts, random_allele=random_allele,
                         folder_out=folder_out, prefix_out_data=prefix_out, 
                         excluded=exclude_pops, diploid_ref=diploid_ref)
    
    ### Set the paths to ref & target
    hmm.p_obj.set_params(h5_path1000g = h5_path1000g, path_targets = path_targets, 
                         meta_path_ref = meta_path_ref, n_ref=n_ref)
    hmm.load_data(iid=iid, ch=ch)  # Load the actual Data
    hmm.load_secondary_objects(c=c)
    
    ### Set the Parameters
    hmm.e_obj.set_params(e_rate = e_rate, e_rate_ref = e_rate_ref)
    hmm.t_obj.set_params(roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump)
    return hmm.compute_tot_neg_likelihood(c)

def hapsb_multiChunk(c, chunks, iid, path_targets_prefix, h5_path1000g, meta_path_ref,
                folder_out, conPop=["CEU"], roh_in=1, roh_out=0, roh_jump=300, e_rate=0.01, e_rate_ref=1e-3,
                processes=1, save=False, save_fp=False, n_ref=2504, diploid_ref=True, 
                exclude_pops=[], e_model="readcount_contam", p_model="SardHDF5", 
                readcounts=True, random_allele=False, prefix_out="", logfile=False):
    # chunks is a dictionary: chrom -> (start of ROH, end of ROH)
    print(f'contamination rate: {c}')
    tot_neg_loglik = 0
    if processes == 1:
        print(f'running using single process...')
        for ch, (start, end) in chunks.items():
            path_targets = path_targets_prefix + f"{iid}.chr{ch}.hdf5"
            tot_neg_loglik += hapsb_chunk_negloglik(iid, ch, start, end, path_targets, h5_path1000g, meta_path_ref,
                folder_out, c, conPop=conPop, roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump, e_rate=e_rate, e_rate_ref=e_rate_ref,
                save=save, save_fp=save_fp, n_ref=n_ref, diploid_ref=diploid_ref, 
                exclude_pops=exclude_pops, e_model=e_model, p_model=p_model, 
                readcounts=readcounts, random_allele=random_allele, prefix_out=prefix_out, logfile=logfile)
    else:
        print(f'running using {processes} processes...')
        prms = [ [iid, ch, start, end, path_targets_prefix + f"{iid}.chr{ch}.hdf5", 
                h5_path1000g, meta_path_ref, folder_out, 
                c, conPop, roh_in, roh_out, roh_jump, e_rate, e_rate_ref,
                save, save_fp, n_ref, diploid_ref, exclude_pops, e_model, p_model, 
                readcounts, random_allele, prefix_out, logfile] 
                for ch, (start, end) in chunks.items()]
        results = multi_run(hapsb_chunk_negloglik, prms, processes = processes)
        print(f'results is: {results}')
        if isinstance(results, float):
            tot_neg_loglik = results
        else:
            tot_neg_loglik += sum(results)
    return tot_neg_loglik

def hapsb_femaleROHcontam(iid, roh_list, path_targets_prefix, h5_path1000g, meta_path_ref,
                folder_out, init_c=0.025, trim=0.0025, minLen=0.045, conPop=["CEU"], roh_in=1, roh_out=0, roh_jump=300, e_rate=0.01, e_rate_ref=1e-3,
                processes=1, save=False, save_fp=False, n_ref=2504, diploid_ref=True, 
                exclude_pops=[], e_model="readcount_contam", p_model="SardHDF5", 
                readcounts=True, random_allele=False, prefix_out="", logfile=False):
    chunks = {}
    with open(roh_list) as f:
        f.readline()
        line = f.readline()
        while line:
            _, _, StartM, EndM, _, lengthM, _, ch, _, _ = line.strip().split(',')
            StartM, EndM, lengthM = float(StartM), float(EndM), float(lengthM)
            if lengthM >= minLen:
                chunks[ch] = (StartM + trim, EndM - trim)
            line = f.readline()
            
    if len(chunks) > 0:
        print(f'a total of {len(chunks)} ROH blocks found.')
        if not path_targets_prefix.endswith('/'):
            path_targets_prefix += "/"
        kargs = (chunks, iid, path_targets_prefix, h5_path1000g, meta_path_ref, folder_out,
                conPop, roh_in, roh_out, roh_jump, e_rate, e_rate_ref, processes, save, save_fp, n_ref, diploid_ref, 
                exclude_pops, e_model, p_model, readcounts, random_allele, prefix_out, logfile)
        res = minimize(hapsb_multiChunk, init_c, args=kargs, method='L-BFGS-B', bounds=[(0, 0.5)])
        if res.success:
            Hfun = ndt.Hessian(hapsb_multiChunk, step=1e-4, full_output=True)
            h, info = Hfun(res.x[0], *kargs)
            h = h[0][0]
            se = math.sqrt(1/(h))
            return res.x[0], se
        else:
            print('L-BFGS-B does not converge...')
            print('Cannot estimate contamination rate from the given data, maybe coverage is too low')
            sys.exit()
    else:
        print(f'not enough ROH blocks found to estimate contamination...')
        sys.exit()


def hapsb_chrom(iid, ch=3, save=True, save_fp=False, n_ref=2504, diploid_ref=True, exclude_pops=[], 
                e_model="EigenstratPacked", p_model="MosaicHDF5", readcounts=True, random_allele=True,
                post_model="Standard", path_targets = "./Data/SA_1240kHDF5/IPK12.h5",
                h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/all1240/chr", 
                meta_path_ref = "./Data/1000Genomes/Individuals/meta_df_all.csv",
                folder_out="./Empirical/Eigenstrat/Reichall/test/", prefix_out="",
                c=0.0, conPop=["CEU"], roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.0,
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
    diploid_ref: Use both haplotypes of reference panel [bool]
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
    parameters = locals() # Gets dictionary of all local variables at this point
    
    ### Create Folder if needed, and pipe output if wanted
    _ = prepare_path(folder_out, iid, ch, prefix_out, logfile=logfile) # Set the logfile
    hmm = HMM_Analyze(cython=3, p_model=p_model, e_model=e_model, post_model=post_model,
                      manual_load=True, save=save, save_fp=save_fp)

    ### Load and prepare the pre-processing Model
    hmm.load_preprocessing_model(conPop)              # Load the preprocessing Model
    hmm.p_obj.set_params(readcounts = readcounts, random_allele=random_allele,
                         folder_out=folder_out, prefix_out_data=prefix_out, 
                         excluded=exclude_pops, diploid_ref=diploid_ref)
    
    ### Set the paths to ref & target
    hmm.p_obj.set_params(h5_path1000g = h5_path1000g, path_targets = path_targets, 
                         meta_path_ref = meta_path_ref, n_ref=n_ref)
    hmm.load_data(iid=iid, ch=ch)  # Load the actual Data
    hmm.load_secondary_objects(c=c)
    
    ### Print out the Parameters used in run:
    print("\nParameters in hapsb_chrom:")
    print("\n".join("{}\t{}".format(k, v) for k, v in parameters.items()))
    print("\n")

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
              path_targets_prefix="", path_targets = "",
              h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/all1240int8/chr", 
              meta_path_ref = "./Data/1000Genomes/Individuals/meta_df_all.csv",
              folder_out="./Empirical/Eigenstrat/Reichall/test/", prefix_out="",
              e_model="haploid", p_model="Eigenstrat", post_model="Standard",
              processes=1, delete=False, output=True, save=True, save_fp=False, 
              n_ref=2504, diploid_ref=True, exclude_pops=[], readcounts=True, random_allele=True,
              c=0.0, conPop=["CEU"], roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.00, 
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
    if len(path_targets) != 0:
        prms = [[iid, ch, save, save_fp, n_ref, diploid_ref, exclude_pops, e_model, p_model, readcounts, random_allele,
            post_model, path_targets, h5_path1000g, meta_path_ref, folder_out, prefix_out,
            c, conPop, roh_in, roh_out, roh_jump, e_rate, e_rate_ref, max_gap, cutoff_post, roh_min_l, logfile] for ch in chs]
    elif len(path_targets_prefix) != 0:
        prms = [[iid, ch, save, save_fp, n_ref, diploid_ref, exclude_pops, e_model, p_model, readcounts, random_allele,
            post_model, f'{path_targets_prefix}/{iid}.chr{ch}.hdf5', h5_path1000g, meta_path_ref, folder_out, prefix_out,
            c, conPop, roh_in, roh_out, roh_jump, e_rate, e_rate_ref, max_gap, cutoff_post, roh_min_l, logfile] for ch in chs]
    else:
        print(f'You need to at least specify one of path_targets or path_targets_prefix...')
        sys.exit()
    assert(len(prms[0])==28)   # Sanity Check
                            
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
        
        
############################################################################
### Run multiple X Chromosomes in parallel

def hapsb_chromXs(iids=[["I15595","I15970"]], ch=23, processes=1, 
                  path_targets = "/project2/jnovembre/hringbauer/caribbean_roh/data/eigenstrat/v421_CaribIllu1000GancSam_bySite_PAM",
                  h5_path1000g = "/project2/jnovembre/hringbauer/HAPSBURG/Data/1000Genomes/HDF5/1240kHDF5/all1240/chr", 
                  meta_path_ref = "/project2/jnovembre/hringbauer/HAPSBURG/Data/1000Genomes/Individuals/meta_df_all_sex.tsv",
                  folder_out = "/project2/jnovembre/hringbauer/HAPSBURG/Empirical/dumpster/testx/", prefix_out="",
                  e_model="readcount", p_model="EigenstratX", post_model="IBD_X",
                  delete=False, output=True, save=True, save_fp=False, 
                  n_ref=2504, diploid_ref=False, exclude_pops=[], readcounts=True, random_allele=False,
                  roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.00, 
                  cutoff_post = 0.999, max_gap=0, roh_min_l = 0.01, logfile=True, combine=False, 
                  file_result="_roh_full.csv"):
    """Run multiple X chromosome pairs. iids is list of two individuals.
    For parameters see hapsb_chrom. Additional:
    iids: list of pairs of iids to run [list of lists of len 2]
    processes: How many processes to use [int]. Think about sufficient memory.
    """
    ### Sanity Checks
    assert(len(iids)>=1)
    assert(len(iids[0])==2)
    
    ### Prepare the Parameters for each Individual pairs
    prms = [[iid, ch, save, save_fp, n_ref, diploid_ref, exclude_pops, e_model, p_model, readcounts, random_allele,
            post_model, path_targets, h5_path1000g, meta_path_ref, folder_out, prefix_out,
            roh_in, roh_out, roh_jump, e_rate, e_rate_ref, max_gap, cutoff_post, roh_min_l, logfile] for iid in iids]
    assert(len(prms[0])==26)   # Sanity Check
                            
    ### Run the analysis in parallel
    multi_run(hapsb_chrom, prms, processes = processes)  
    
    ### Move results to main folder if needed
    if len(file_result)>0:
        for iid in iids:
            move_X_to_parent_folder(base_path=folder_out, 
                                    iid=iid, delete=delete, ch=ch, 
                                    prefix_out=prefix_out, file_result=file_result)
    
    
    
##################################################################################
##################################################################################
######################## contamination code starts here ##########################

def hapCon_chrom(iid, ch, save=True, save_fp=False, n_ref=2504, diploid_ref=True, exclude_pops=[], conPop=[], 
                e_model="readcount_contam", p_model="SardHDF5", readcounts=True, random_allele=False,
                post_model="Standard", path_targets = "./Data/SA_1240kHDF5/IPK12.h5",
                h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/all1240/chr", 
                meta_path_ref = "./Data/1000Genomes/Individuals/meta_df_all.csv",
                folder_out="./Empirical/Eigenstrat/Reichall/test/", prefix_out="",
                c=np.arange(0, 0.2, 0.005), roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.0,
                max_gap=0, cutoff_post = 0.999, roh_min_l = 0.01, logfile=True):
    """Run HapCon analysis for one chromosome on hdf5 data
    Wrapper for HMM Class.
    iid: IID of the Target Individual, as found in the given hdf5 file [str]
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
    diploid_ref: Use both haplotypes of reference panel [bool]
    exclude_pops: Which populations to exclude from reference [list of str]
    conPop: use which population in the ref panel as the contaminating pop. If empty list, then use all samples in the ref panel to cauclate allele freq.
    readcounts: Whether to load readcount data [bool]
    random_allele: Whether to pick a random of the two target alleles per locus [bool]
    c: contamination rate
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
    
    RETURN: loglikelihood of the data under the model
    
    """
    parameters = locals() # Gets dictionary of all local variables at this point
    
    ### Create Folder if needed, and pipe output if wanted
    _ = prepare_path(folder_out, iid, ch, prefix_out, logfile=logfile) # Set the logfile
    hmm = HMM_Analyze(cython=3, p_model=p_model, e_model=e_model, post_model=post_model,
                      manual_load=True, save=save, save_fp=save_fp)

    ### Load and prepare the pre-processing Model
    hmm.load_preprocessing_model(conPop)              # Load the preprocessing Model
    hmm.p_obj.set_params(readcounts = readcounts, random_allele=random_allele,
                         folder_out=folder_out, prefix_out_data=prefix_out, 
                         excluded=exclude_pops, diploid_ref=diploid_ref)
    
    ### Set the paths to ref & target
    hmm.p_obj.set_params(h5_path1000g = h5_path1000g, path_targets = path_targets, 
                         meta_path_ref = meta_path_ref, n_ref=n_ref)
    hmm.load_data(iid=iid, ch=ch)  # Load the actual Data
    hmm.load_secondary_objects()
    
    ### Print out the Parameters used in run:
    print("\nParameters in hapCon_chrom:")
    print("\n".join("{}\t{}".format(k, v) for k, v in parameters.items()))
    print("\n")

    ### Set the Parameters
    hmm.e_obj.set_params(e_rate = e_rate, e_rate_ref = e_rate_ref)
    hmm.t_obj.set_params(roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump)
    hmm.post_obj.set_params(max_gap=max_gap, cutoff_post=cutoff_post, roh_min_l = roh_min_l)


    lls, con_mle, lower, upper = hmm.optimize_ll_contamination(c)
    return lls, con_mle, lower, upper

    
    # calculate the posterior. Should set full=True later as there is no need to store posterior in a file for contamination purpose
    # for debugging purpose now, leave full=False as is defaulted
    #*_, tot_ll = hmm.calc_posterior(save=save, full=True) # Calculate the Posterior.
    #hmm.calc_posterior(save=save)
    # Do the Post-Processing. Just here for sanity check of called ROH region. 
    # Not needed for estimating contamination purpose. TO BE REMOVED LATER.
    #hmm.post_processing(save=save)             
    #print(f'hapcon_chr returns result: {tot_ll}')
    #return tot_ll

def hapCon_ind(iid, chs=range(1,23), 
            path_targets='/mnt/archgen/users/yilei/Data/SA_1240KHDF5/marcus2020.h5',
            h5_path1000g='./Data/1000G/1000g1240khdf5/all1240/chr',
            meta_path_ref='./Data/1000G/meta_df_all.csv', 
            folder_out='./tools/hapROH/test_output/contaminate/first_try/', prefix_out="",
            e_model="readcount_contam", p_model="SardHDF5", post_model="Standard",
            processes=1, delete=False, output=True, save=True, save_fp=False, 
            n_ref=2504, diploid_ref=True, exclude_pops=[], conPop=[], readcounts=True, random_allele=False,
            c=0.05, roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.00, 
            cutoff_post = 0.999, max_gap=0, roh_min_l = 0.01, logfile=False):
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
    c: contamination rate
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
    prms = [[iid, ch, save, save_fp, n_ref, diploid_ref, exclude_pops, conPop, e_model, p_model, readcounts, random_allele,
            post_model, path_targets, h5_path1000g, meta_path_ref, folder_out, prefix_out, c, 
            roh_in, roh_out, roh_jump, e_rate, e_rate_ref, max_gap, cutoff_post, roh_min_l, logfile] for ch in chs]
    assert(len(prms[0]) == 28)   # Sanity Check
                            
    ### Run the analysis in parallel
    #results = multi_run(hapCon_chrom, prms, processes = processes)

    # just for now, do this sequentially
    lls = []
    for prm in prms:
        lls.append(hapCon_chrom(*prm))

    tot_lls = sum(lls)
    print(f'total log likelihood: {tot_lls}')
    return tot_lls



def hapCon_chrom_2d(iid, ch, save=True, save_fp=False, n_ref=2504, diploid_ref=True, exclude_pops=[], conPop=[], 
                e_model="readcount_contam", p_model="SardHDF5", readcounts=True, random_allele=False,
                post_model="Standard", path_targets = "./Data/SA_1240kHDF5/IPK12.h5",
                h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/all1240/chr", 
                meta_path_ref = "./Data/1000Genomes/Individuals/meta_df_all.csv",
                folder_out="./Empirical/Eigenstrat/Reichall/test/", prefix_out="",
                c=0.025, roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.0,
                max_gap=0, cutoff_post = 0.999, roh_min_l = 0.01, logfile=True):
    """Run HapCon analysis for one chromosome on hdf5 data
    Wrapper for HMM Class.
    iid: IID of the Target Individual, as found in the given hdf5 file [str]
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
    diploid_ref: Use both haplotypes of reference panel [bool]
    exclude_pops: Which populations to exclude from reference [list of str]
    conPop: use which population in the ref panel as the contaminating pop. If empty list, then use all samples in the ref panel to cauclate allele freq.
    readcounts: Whether to load readcount data [bool]
    random_allele: Whether to pick a random of the two target alleles per locus [bool]
    c: contamination rate
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
    
    RETURN: loglikelihood of the data under the model
    
    """
    parameters = locals() # Gets dictionary of all local variables at this point
    
    ### Create Folder if needed, and pipe output if wanted
    _ = prepare_path(folder_out, iid, ch, prefix_out, logfile=logfile) # Set the logfile
    hmm = HMM_Analyze(cython=3, p_model=p_model, e_model=e_model, post_model=post_model,
                      manual_load=True, save=save, save_fp=save_fp)

    ### Load and prepare the pre-processing Model
    hmm.load_preprocessing_model(conPop)              # Load the preprocessing Model
    hmm.p_obj.set_params(readcounts = readcounts, random_allele=random_allele,
                         folder_out=folder_out, prefix_out_data=prefix_out, 
                         excluded=exclude_pops, diploid_ref=diploid_ref)
    
    ### Set the paths to ref & target
    hmm.p_obj.set_params(h5_path1000g = h5_path1000g, path_targets = path_targets, 
                         meta_path_ref = meta_path_ref, n_ref=n_ref)
    hmm.load_data(iid=iid, ch=ch)  # Load the actual Data
    hmm.load_secondary_objects()
    
    ### Print out the Parameters used in run:
    print("\nParameters in hapCon_chrom:")
    print("\n".join("{}\t{}".format(k, v) for k, v in parameters.items()))
    print("\n")

    ### Set the Parameters
    hmm.e_obj.set_params(e_rate = e_rate, e_rate_ref = e_rate_ref)
    hmm.t_obj.set_params(roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump)
    hmm.post_obj.set_params(max_gap=max_gap, cutoff_post=cutoff_post, roh_min_l = roh_min_l)


    mle, se = hmm.optimize_ll_contamANDerr(c, e_rate)
    return mle, se


def hapCon_chrom_BFGS(iid, ch, save=True, save_fp=False, n_ref=2504, diploid_ref=True, exclude_pops=[], conPop=[], 
    e_model="readcount_contam", p_model="SardHDF5", readcounts=True, random_allele=False,
    post_model="Standard", path_targets = "./Data/SA_1240kHDF5/IPK12.h5",
    h5_path1000g = "./Data/1000Genomes/HDF5/1240kHDF5/all1240/chr", 
    meta_path_ref = "./Data/1000Genomes/Individuals/meta_df_all.csv",
    folder_out="./Empirical/Eigenstrat/Reichall/test/", prefix_out="",
    c=0.025, roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.0,
    max_gap=0, cutoff_post = 0.999, roh_min_l = 0.01, logfile=True):

    parameters = locals() # Gets dictionary of all local variables at this point
    
    ### Create Folder if needed, and pipe output if wanted
    _ = prepare_path(folder_out, iid, ch, prefix_out, logfile=logfile) # Set the logfile
    hmm = HMM_Analyze(cython=3, p_model=p_model, e_model=e_model, post_model=post_model,
                      manual_load=True, save=save, save_fp=save_fp)

    ### Load and prepare the pre-processing Model
    hmm.load_preprocessing_model(conPop)              # Load the preprocessing Model
    hmm.p_obj.set_params(readcounts = readcounts, random_allele=random_allele,
                         folder_out=folder_out, prefix_out_data=prefix_out, 
                         excluded=exclude_pops, diploid_ref=diploid_ref)
    
    ### Set the paths to ref & target
    hmm.p_obj.set_params(h5_path1000g = h5_path1000g, path_targets = path_targets, 
                         meta_path_ref = meta_path_ref, n_ref=n_ref)
    hmm.load_data(iid=iid, ch=ch)  # Load the actual Data
    print(hmm.pCon)
    hmm.load_secondary_objects()
    
    ### Print out the Parameters used in run:
    print("\nParameters in hapCon_chrom_BFGS:")
    print("\n".join("{}\t{}".format(k, v) for k, v in parameters.items()))
    print("\n")

    ### Set the Parameters
    hmm.e_obj.set_params(e_rate = e_rate, e_rate_ref = e_rate_ref)
    hmm.t_obj.set_params(roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump)
    hmm.post_obj.set_params(max_gap=max_gap, cutoff_post=cutoff_post, roh_min_l = roh_min_l)


    con_mle, lower, upper = hmm.optimize_ll_contamination_BFGS(c)
    return con_mle, lower, upper



##################################################################################
##################################################################################
################################### END ##########################################


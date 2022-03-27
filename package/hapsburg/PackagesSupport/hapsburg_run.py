"""
Function to run HAPSBURG on a full individual and full reference Dataset,
with all relevant keywords.
Function for running on single chromsome, with all relevant keywords.
@Author: Harald Ringbauer, 2019
"""

import math
import numpy as np
import sys
import os
import time
from scipy.optimize import minimize
import numdifftools as ndt
from scipy.optimize import newton
import shutil



from hapsburg.hmm_inference import HMM_Analyze   # The HMM core object
from hapsburg.PackagesSupport.parallel_runs.helper_functions import prepare_path, multi_run, combine_individual_data, move_X_to_parent_folder
from hapsburg.PackagesSupport.loadEigenstrat.saveHDF5 import mpileup2hdf5, bam2hdf5, bamTable2hdf5, mpileup2hdf5_damageAware

# def hapsb_chunk_negloglik(iid, ch, start, end, path_targets, h5_path1000g, meta_path_ref,
#                 folder_out, c, conPop=["CEU"], roh_in=1, roh_out=0, roh_jump=300, e_rate=0.01, e_rate_ref=1e-3,
#                 save=False, save_fp=False, n_ref=2504, diploid_ref=True, 
#                 exclude_pops=[], e_model="readcount_contam", p_model="SardHDF5", 
#                 readcounts=True, random_allele=False, prefix_out="", logfile=False):
#     parameters = locals() # Gets dictionary of all local variables at this point
    
#     ### Create Folder if needed, and pipe output if wanted
#     _ = prepare_path(folder_out, iid, ch, prefix_out, logfile=logfile) # Set the logfile
#     hmm = HMM_Analyze(cython=3, p_model=p_model, e_model=e_model, output=False,
#                       manual_load=True, save=save, save_fp=save_fp, start=start, end=end)

#     ### Load and prepare the pre-processing Model
#     hmm.load_preprocessing_model(conPop)              # Load the preprocessing Model
#     hmm.p_obj.set_params(readcounts = readcounts, random_allele=random_allele,
#                          folder_out=folder_out, prefix_out_data=prefix_out, 
#                          excluded=exclude_pops, diploid_ref=diploid_ref)
    
#     ### Set the paths to ref & target
#     hmm.p_obj.set_params(h5_path1000g = h5_path1000g, path_targets = path_targets, 
#                          meta_path_ref = meta_path_ref, n_ref=n_ref)
#     hmm.load_data(iid=iid, ch=ch)  # Load the actual Data
#     hmm.load_secondary_objects(c=c)
    
#     ### Set the Parameters
#     hmm.e_obj.set_params(e_rate = e_rate, e_rate_ref = e_rate_ref)
#     hmm.t_obj.set_params(roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump)
#     return hmm.compute_tot_neg_likelihood(c)

def prepare_path_general(basepath, iid, prefix, suffix, logfile):
    if not os.path.exists(basepath):
        os.makedirs(basepath)
    ### Activate LOG FILE output if given
    if logfile == True:
        if prefix:
            path_log = os.path.join(basepath, f"{iid}.{prefix}.{suffix}.log")
        else:
            path_log = os.path.join(basepath, f"{iid}.{suffix}.log")
        print(f"Set Output Log path: {path_log}")
        sys.stdout = open(path_log, 'w')

def preload(iid, ch, start, end, path_targets, h5_path1000g, meta_path_ref,
                folder_out, conPop=["CEU"], roh_jump=300, e_rate=0.01, e_rate_ref=1e-3,
                n_ref=2504, exclude_pops=[]):
    # reference panel needed only to be loaded once for each optimization phase
    # so we preload here to reduce run time
        
    hmm = HMM_Analyze(cython=3, p_model="SardHDF5", e_model="readcount_contam", output=False,
                      manual_load=True, save=False, save_fp=False, start=start, end=end)

    ### Load and prepare the pre-processing Model
    hmm.load_preprocessing_model(conPop)              # Load the preprocessing Model
    hmm.p_obj.set_params(readcounts=True, random_allele=False,
                         folder_out=folder_out, prefix_out_data="", 
                         excluded=exclude_pops, diploid_ref=True)
    
    ### Set the paths to ref & target
    hmm.p_obj.set_params(h5_path1000g = h5_path1000g, path_targets = path_targets, 
                         meta_path_ref = meta_path_ref, n_ref=n_ref)
    hmm.load_data(iid=iid, ch=ch)  # Load the actual Data
    hmm.load_secondary_objects()
    
    ### Set the Parameters
    hmm.e_obj.set_params(e_rate = e_rate, e_rate_ref = e_rate_ref)
    hmm.t_obj.set_params(roh_out=0.0, roh_jump=roh_jump)
    return hmm
    

# def hapsb_multiChunk(c, chunks, iid, path_targets_prefix, h5_path1000g, meta_path_ref,
#                 folder_out, conPop=["CEU"], roh_in=1, roh_out=0, roh_jump=300, e_rate=0.01, e_rate_ref=1e-3,
#                 processes=1, save=False, save_fp=False, n_ref=2504, diploid_ref=True, 
#                 exclude_pops=[], e_model="readcount_contam", p_model="SardHDF5", 
#                 readcounts=True, random_allele=False, prefix_out="", logfile=False):
#     # chunks is a dictionary: chrom -> (start of ROH, end of ROH)
#     t1 = time.time()
#     tot_neg_loglik = 0
#     if processes == 1:
#         # print(f'running using single process...')
#         for ch, start, end in chunks:
#             path_targets = path_targets_prefix + f"{iid}.chr{ch}.hdf5"
#             tot_neg_loglik += hapsb_chunk_negloglik(iid, ch, start, end, path_targets, h5_path1000g, meta_path_ref,
#                 folder_out, c, conPop=conPop, roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump, e_rate=e_rate, e_rate_ref=e_rate_ref,
#                 save=save, save_fp=save_fp, n_ref=n_ref, diploid_ref=diploid_ref, 
#                 exclude_pops=exclude_pops, e_model=e_model, p_model=p_model, 
#                 readcounts=readcounts, random_allele=random_allele, prefix_out=prefix_out, logfile=logfile)
#     else:
#         # print(f'running using {processes} processes...')
#         prms = [ [iid, ch, start, end, path_targets_prefix + f"{iid}.chr{ch}.hdf5", 
#                 h5_path1000g, meta_path_ref, folder_out, 
#                 c, conPop, roh_in, roh_out, roh_jump, e_rate, e_rate_ref,
#                 save, save_fp, n_ref, diploid_ref, exclude_pops, e_model, p_model, 
#                 readcounts, random_allele, prefix_out, logfile] 
#                 for ch, start, end in chunks]
#         results = multi_run(hapsb_chunk_negloglik, prms, processes = processes)
#         # print(f'results is: {results}')
#         if isinstance(results, float):
#             tot_neg_loglik = results
#         else:
#             tot_neg_loglik += sum(results)
#     print(f'hapsb_multichunk takes {time.time()-t1}')
#     return tot_neg_loglik

def hapsb_chunk_negloglik_preload(hmm, c):
    ret = hmm.compute_tot_neg_likelihood(c)
    return ret

def hapsb_multiChunk_preload(c, hmms, processes=1):
    t1 = time.time()
    tot_neg_loglik = 0
    if processes == 1:
        for hmm in hmms:
            tot_neg_loglik += hmm.compute_tot_neg_likelihood(c)
    else:
        prms = [[hmm, c] for hmm in hmms ]
        results = multi_run(hapsb_chunk_negloglik_preload, prms, processes=processes)
        if isinstance(results, float):
            tot_neg_loglik = results
        else:
            tot_neg_loglik += sum(results)
    print(f"hapsb_multiChunk_preload takes {time.time()-t1}")
    return tot_neg_loglik



# def hapsb_femaleROHcontam(iid, roh_list, path_targets_prefix, h5_path1000g, meta_path_ref,
#                 folder_out, init_c=0.025, trim=0.005, minLen=0.05, conPop=["CEU"], roh_in=1, roh_out=0, roh_jump=300, e_rate=0.01, e_rate_ref=1e-3,
#                 processes=1, save=False, save_fp=False, n_ref=2504, diploid_ref=True, 
#                 exclude_pops=[], e_model="readcount_contam", p_model="SardHDF5", 
#                 readcounts=True, random_allele=False, prefix_out="", logfile=False):
#     chunks = []
#     with open(roh_list) as f:
#         f.readline()
#         line = f.readline()
#         while line:
#             _, _, StartM, EndM, _, lengthM, _, ch, _, _ = line.strip().split(',')
#             StartM, EndM, lengthM = float(StartM), float(EndM), float(lengthM)
#             if lengthM >= minLen:
#                 chunks.append((ch, StartM + trim, EndM - trim))
#                 print(f'chr{ch}\t{round(StartM, 6)}\t{round(EndM, 6)}')
#             line = f.readline()
            
#     if len(chunks) > 0:
#         sumROH = 0
#         for _, start, end in chunks:
#             sumROH += end - start
#         print(f'a total of {len(chunks)} ROH blocks passing filtering threshold found, total length after trimming: {sumROH}M.')
#         if not path_targets_prefix.endswith('/'):
#             path_targets_prefix += "/"
#         kargs = (chunks, iid, path_targets_prefix, h5_path1000g, meta_path_ref, folder_out,
#                 conPop, roh_in, roh_out, roh_jump, e_rate, e_rate_ref, processes, save, save_fp, n_ref, diploid_ref, 
#                 exclude_pops, e_model, p_model, readcounts, random_allele, prefix_out, logfile)
#         res = minimize(hapsb_multiChunk, init_c, args=kargs, method='L-BFGS-B', bounds=[(0, 0.5)])
#         if not res.success:
#             print('L-BFGS-B does not converge. Printing its result log for diagnostic purpose.')
#             print(res)
#             print('please take the final estimate with caution.')
#         Hfun = ndt.Hessian(hapsb_multiChunk, step=1e-4, full_output=True)
#         h, info = Hfun(res.x[0], *kargs)
#         h = h[0][0]
#         se = math.sqrt(1/(h))
#         return res.x[0], se
#     else:
#         print(f'not enough ROH blocks found to estimate contamination...')
#         sys.exit()

def hapsb_femaleROHcontam_preload(iid, roh_list, h5_path1000g, meta_path_ref,
                hdf5_path=None, folder_out=None, 
                init_c=0.025, trim=0.005, minLen=0.05, conPop=["CEU"], roh_jump=300, 
                e_rate=1e-2, e_rate_ref=1e-3, processes=1, n_ref=2504, 
                exclude_pops=["AFR"], prefix=None, logfile=False, cleanup=False):
    """
    Estimating autosomal contamination rate from a list of ROH blocks. Need at least one ROH for inference.

    Parameters
    ----------
    iid: str
        IID of the sample. We assume that the mpileup file has the format $iid.chr[1-22].mpileup.
    roh_list: str
        Path to a file containing a list of ROH blocks. This file should have the same format as the output of hapROH.
    h5_path1000g: str
        Path to the reference panel.
    meta_path_ref: str
        Path to the metadata of reference panel.
    hdf5_path: str
        Directory of hdf5 files. One file for each autosome.
    folder_out: str
        Directory in which you want the output to reside. If not given, all output files will be in the parent directory of mpileup_path.
    init_c: float
        Initial value for the BFGS search.
    trim: float
        Trim both ends of inferred ROH blocks (in Morgan).
    minLen: float
        Minimum length of ROH blocks to use in estimating contamination (in Morgan).
    conPop: list of str
        Contaminant Ancestry. Must correspond to names in the super_pop or pop column in the 1000G metadata file.
    roh_jump: float
        Copying jump rate.
    e_rate: float
        If mpileup_path is provided, the sequencing error rate will be estimated from flanking sites, so this parameter has no effect.
        If hdf5_path is provided, then this parameter is the error rate used for inference. So one should obtain sequencing error estimate from external source if you use the hdf5_path option.
    e_rate_ref: float
        Haplotype copying error rate.
    processes: int
        Number of processes to use.
    n_ref: int
        Number of samples in the reference panel.
    exclude_pops: list of str
        A list of populations to exclude from the reference panel.
    prefix: str
        Prefix of the output and log file. The output will follow $iid.$prefix.hapCON_ROH.txt. And the log file will follow $iid.$prefix.hapCON_ROH.log.
    logfile: bool
        Whether to produce a log file.
    cleanup: bool
        Whether to delete intermediary HDF5 files generated during this function run.

    Returns
    ---------
    conMLE: float
        MLE estimate for contamination.
    se: float
        Standard error of the estimated contamination rate.
     
    """
    
    
    # should be the same as hapsb_femaleROHcontam, but a faster implementation
    if not folder_out:
        folder_out = os.path.dirname(os.path.abspath(hdf5_path))

    if prefix:
        fileName = f'{iid}.{prefix}.hapCON_ROH.txt'
    else:
        fileName = f'{iid}.hapCON_ROH.txt'

    prepare_path_general(folder_out, iid, prefix, "hapCON_ROH", logfile)
    chunks = []
    with open(roh_list) as f:
        f.readline()
        line = f.readline()
        while line:
            _, _, StartM, EndM, _, lengthM, _, ch, _, _ = line.strip().split(',')
            StartM, EndM, lengthM = float(StartM), float(EndM), float(lengthM)
            if lengthM >= minLen:
                chunks.append((ch, StartM + trim, EndM - trim))
                print(f'chr{ch}\t{round(StartM, 6)}\t{round(EndM, 6)}')
            line = f.readline()
            
    if len(chunks) > 0:
        sumROH = 0
        for _, start, end in chunks:
            sumROH += end - start
        print(f'a total of {len(chunks)} ROH blocks passing filtering threshold found, total length after trimming: {100*sumROH:.3f}cM.')


        # preload hmm models
        t1 = time.time()
        prms = [[iid, ch, start, end, hdf5_path+f"/{iid}.chr{ch}.hdf5", h5_path1000g, meta_path_ref, \
                    folder_out, conPop, roh_jump, e_rate, e_rate_ref, n_ref, exclude_pops] for ch, start, end in chunks]
        hmms = multi_run(preload, prms, processes=processes)
        assert(len(hmms) == len(chunks))
        # hmms = []
        # for ch, start, end in chunks:
        #     path_targets = hdf5_path + "/" +f"{iid}.chr{ch}.hdf5"
        #     hmm = preload(iid, ch, start, end, path_targets, h5_path1000g, meta_path_ref,
        #         folder_out, conPop=conPop, roh_jump=roh_jump, 
        #         e_rate=e_rate, e_rate_ref=e_rate_ref, n_ref=n_ref, 
        #         exclude_pops=exclude_pops)
        #     hmms.append(hmm)
        print(f'{len(chunks)} hmm models loaded, takes {round(time.time()-t1, 3)}s')

        # the actual optimization part
        kargs = (hmms, processes)
        res = minimize(hapsb_multiChunk_preload, init_c, args=kargs, method='L-BFGS-B', bounds=[(0, 0.5)])
        if not res.success:
            print('L-BFGS-B does not converge. Printing its result log for diagnostic purpose.')
            print(res)
            print('please treat the final estimate with caution.')
        Hfun = ndt.Hessian(hapsb_multiChunk_preload, step=1e-4, full_output=True)
        try:
            x = res.x[0]
            h, info = Hfun(x, *kargs)
            h = h[0][0]
            if h < 0:
                print('WARNING: Cannot estimate standard error because the likelihood curve is concave up...')
                se = np.nan
            else:
                if x > 0:
                    se = math.sqrt(1/(h))
                else:
                    # hessian does not work well at the boundary, use a different approach
                    print(f'use quadracitc interpolation to obtain likelihood confidence interval...')
                    step = 1e-6
                    grad = (hapsb_multiChunk_preload(step, *kargs) - hapsb_multiChunk_preload(0, *kargs))/step
                    assert(grad > 0)
                    findroot = lambda x, x0, grad, hess: hess*(x-x0)**2/2.0 + (x-x0)*grad - 1.92
                    findroot_prime = lambda x, x0, grad, hess: (x-x0)*hess + grad
                    res = newton(findroot, x, fprime=findroot_prime, args=(x, grad, h))
                    se = res/1.96
            # with open(f'{folder_out}/{fileName}', 'w') as f:
            #     f.write(f'Method1: Fixing genotyping error rate at {e_rate}\n')
            #     f.write(f'\tROH blocks obtained from: {roh_list}\n')
            #     f.write(f'\tNumber of ROH blocks found: {len(chunks)}\n')
            #     f.write(f'\tTotal length of ROH after trimming: {round(100*sumROH,3)}cM\n')
            #     f.write(f'\tMLE for contamination using BFGS: {round(x, 6)} ({round(x-1.96*se, 6)} - {round(x+1.96*se, 6)})\n')
            return x, se, chunks, sumROH
        except AssertionError:
            print(f'cannot estimate the Hessian of the loglikelihood around {res.x}')
            se = np.nan
        finally:
            if cleanup:
                shutil.rmtree(hdf5_path)
                print(f'deleted intermediary hdf5 files at {hdf5_path}')
            return x, se, chunks, sumROH
    else:
        print(f'not enough ROH blocks found to estimate contamination...')
        with open(f'{folder_out}/{fileName}', 'w') as f:
            f.write('not enough ROH blocks found to estimate contamination...')
        sys.exit()


def hapsb_chrom(iid, ch=3, save=True, save_fp=False, n_ref=2504, diploid_ref=True, exclude_pops=[], 
                e_model="EigenstratPacked", p_model="MosaicHDF5", readcounts=True, random_allele=True,
                post_model="Standard", path_targets=None,
                h5_path1000g=None, meta_path_ref=None, folder_out=None, prefix_out="",
                c=0.0, conPop=["CEU"], roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.0,
                max_gap=0.005, roh_min_l_initial = 0.02, roh_min_l_final = 0.05,
                min_len1 = 0.02, min_len2 = 0.04, cutoff_post = 0.999, logfile=True):
    """Run Hapsburg analysis for one chromosome on eigenstrat data
    Wrapper for HMM Class.

    Parameters
    ----------
    iid: str
        IID of the Target Individual, as found in Eigenstrat.
    ch: int
        Which chromosomes to call ROH.
    path_targets: str
        Path of the target files. You need only specify one of path_targets_prefix or path_targets.
    h5_path1000g: str
        Path of the reference genotypes
    meta_path_ref: str 
        Path of the meta file for the references
    folder_out: str
        Path of the basis folder for output
    prefix_out: str
        Path to insert in output string, e.g. test/ [str]
    e_model: str
        Emission model to use, should be one of haploid/diploid_gt/readcount
    p_model: str
        Preprocessing model to use, should be one of EigenstratPacked/EigenstratUnpacked/MosaicHDF5
    post_model: str
        Model to post-process the data, should be one of Standard/MMR (experimental)
    save: bool
        Whether to save the inferred ROH
    save_fp: bool
        Whether to save the full posterior matrix
    n_ref: int
        Number of (diploid) reference Individuals to use
    diploid: bool
        Whether the reference panel is diploid or not (e.g., for autosome, True and for male X chromosome, False). 
    exclude_pops: list of str
        Which populations to exclude from reference
    readcounts: bool
        Whether to load readcount data
    random_allele: bool
        Whether to pick a random of the two target alleles per locus
    c: float
        Contamination rate. This is only applicable if the emission model is readcount_contam.
    conPop: list of str
        Ancestry of contamination source. Only applicable if the emission model is readcount_contam.
    roh_in: float
        Parater to jump into ROH state (per Morgan)
    roh_out: float
        Parameter to jump out of ROH state (per Morgan)
    roh_jump: float
        Parameter to jump (per Morgan)
    e_rate: float
        Sequencing error rate.
    e_rate_ref: float
        Haplotype miscopying rate.
    cutoff_post: float
        Posterior cutoff for ROH calling
    max_gap: float
        Maximum gap to merge two adjacent short ROH blocks (in Morgan)
    roh_min_l_initial: float
        Minimum length of ROH blocks to use before merging adjacent ones (in Morgan)
    roh_min_l_final: float
        Minimum length of ROH blcoks to output after merging (in Morgan)
    min_len1: float
        Minimum length of the shorter candidate block in two adjacent blocks that can be merged (in Morgan)
    min_len2: float
        Minimum length of the longer candidate block in two adjacent blocks that can be merged (in Morgan)
    logfile: bool
        Whether to use logfile
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
    hmm.load_secondary_objects(c=c)
    
    ### Print out the Parameters used in run:
    print("\nParameters in hapsb_chrom:")
    print("\n".join("{}\t{}".format(k, v) for k, v in parameters.items()))
    print("\n")

    ### Set the Parameters
    hmm.e_obj.set_params(e_rate = e_rate, e_rate_ref = e_rate_ref)
    hmm.t_obj.set_params(roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump)
    hmm.post_obj.set_params(max_gap=max_gap, cutoff_post=cutoff_post,\
                roh_min_l_initial=roh_min_l_initial, roh_min_l_final=roh_min_l_final,\
                min_len1=min_len1, min_len2=min_len2)
    
    ### hmm.calc_viterbi_path(save=save)           # Calculate the Viterbi Path.
    hmm.calc_posterior(save=save)              # Calculate the Posterior.
    hmm.post_processing(save=save)             # Do the Post-Processing.


#########################################################
### Run Hapsburg for one Individual (wrap for Chr.)

def hapsb_ind(iid, chs=range(1,23),
              path_targets_prefix="", path_targets="",
              h5_path1000g=None, meta_path_ref=None, folder_out=None, prefix_out="",
              e_model="haploid", p_model="Eigenstrat", post_model="Standard",
              processes=1, delete=False, output=True, save=True, save_fp=False, 
              n_ref=2504, diploid_ref=True, exclude_pops=[], readcounts=True, random_allele=True,
              c=0.0, conPop=["CEU"], roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.00, 
              cutoff_post = 0.999, max_gap=0.005, roh_min_l_initial = 0.02, roh_min_l_final = 0.05,
                min_len1 = 0.02, min_len2 = 0.04, logfile=True, combine=True, file_result="_roh_full.csv"):
    """Analyze a full single individual in a parallelized fashion. Run multiple chromosome analyses in parallel.
    Then brings together the result ROH tables from each chromosome into one genome-wide summary ROH table.
    This function wraps hapsb_chrom. The default Parameters are finetuned for pseudo-haploid 1240k aDNA data.

    Parameters
    ----------
    iid: str
        IID of the Target Individual, as found in Eigenstrat.
    chs: list
        Which set of chromosomes to call ROH.
    path_targets_prefix:
        A directory containing a hdf5 file for each chromosome. The file name should follow $iid.chr$ch.hdf5.
    path_targets: str
        Path of the target files. You need only specify one of path_targets_prefix or path_targets.
    h5_path1000g: str
        Path of the reference genotypes
    meta_path_ref: str 
        Path of the meta file for the references
    folder_out: str
        Path of the basis folder for output
    prefix_out: str
        Path to insert in output string, e.g. test/ [str]
    e_model: str
        Emission model to use, should be one of haploid/diploid_gt/readcount
    p_model: str
        Preprocessing model to use, should be one of EigenstratPacked/EigenstratUnpacked/MosaicHDF5
    post_model: str
        Model to post-process the data, should be one of Standard/MMR (experimental)
    processes: int
        How many Processes to use
    delete: bool
        Whether to delete raw posterior per locus
    output: bool
        Whether to print extensive output
    save: bool
        Whether to save the inferred ROH
    save_fp: bool
        Whether to save the full posterior matrix
    n_ref: int
        Number of (diploid) reference Individuals to use
    diploid: bool
        Whether the reference panel is diploid or not (e.g., for autosome, True and for male X chromosome, False). 
    exclude_pops: list of str
        Which populations to exclude from reference
    readcounts: bool
        Whether to load readcount data
    random_allele: bool
        Whether to pick a random of the two target alleles per locus
    c: float
        Contamination rate. This is only applicable if the emission model is readcount_contam.
    conPop: list of str
        Ancestry of contamination source. Only applicable if the emission model is readcount_contam.
    roh_in: float
        Parater to jump into ROH state (per Morgan)
    roh_out: float
        Parameter to jump out of ROH state (per Morgan)
    roh_jump: float
        Parameter to jump (per Morgan)
    e_rate: float
        Sequencing error rate.
    e_rate_ref: float
        Haplotype miscopying rate.
    cutoff_post: float
        Posterior cutoff for ROH calling
    max_gap: float
        Maximum gap to merge two adjacent short ROH blocks (in Morgan)
    roh_min_l_initial: float
        Minimum length of ROH blocks to use before merging adjacent ones (in Morgan)
    roh_min_l_final: float
        Minimum length of ROH blcoks to output after merging (in Morgan)
    min_len1: float
        Minimum length of the shorter candidate block in two adjacent blocks that can be merged (in Morgan)
    min_len2: float
        Minimum length of the longer candidate block in two adjacent blocks that can be merged (in Morgan)
    logfile: bool
        Whether to use logfile
    combine: bool 
        Wether to combine output of all chromosomes
    file_result: str
        Appendix to individual results
    """
                            
    if output:
        print(f"Doing Individual {iid}...")
    
    ### Prepare the Parameters for that Indivdiual
    if len(path_targets) != 0:
        prms = [[iid, ch, save, save_fp, n_ref, diploid_ref, exclude_pops, e_model, p_model, readcounts, random_allele,
            post_model, path_targets, h5_path1000g, meta_path_ref, folder_out, prefix_out,
            c, conPop, roh_in, roh_out, roh_jump, e_rate, e_rate_ref, max_gap, roh_min_l_initial, 
            roh_min_l_final, min_len1, min_len2, cutoff_post, logfile] for ch in chs]
    elif len(path_targets_prefix) != 0:
        prms = [[iid, ch, save, save_fp, n_ref, diploid_ref, exclude_pops, e_model, p_model, readcounts, random_allele,
            post_model, f'{path_targets_prefix}/{iid}.chr{ch}.hdf5', h5_path1000g, meta_path_ref, folder_out, prefix_out,
            c, conPop, roh_in, roh_out, roh_jump, e_rate, e_rate_ref, max_gap, roh_min_l_initial, roh_min_l_final, 
            min_len1, min_len2, cutoff_post, logfile] for ch in chs]
    else:
        print(f'You need to at least specify one of path_targets or path_targets_prefix...')
        sys.exit()
    assert(len(prms[0])==31)   # Sanity Check
                            
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

def hapsb_chromXs(iids=[["I15595","I15970"]], processes=1, file_result="", 
                  ch=23, save=True, save_fp=False, n_ref=2504, diploid_ref=True, exclude_pops=[], 
                  e_model="EigenstratPacked", p_model="MosaicHDF5", readcounts=True, random_allele=True,
                  post_model="Standard", path_targets=None,
                  h5_path1000g=None, meta_path_ref=None, folder_out=None, prefix_out="",
                  c=0.0, conPop=["CEU"], roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.0,
                  max_gap=0.005, roh_min_l_initial = 0.02, roh_min_l_final = 0.05,
                  min_len1 = 0.02, min_len2 = 0.04, cutoff_post = 0.999, logfile=True):
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
             c, conPop, roh_in, roh_out, roh_jump, e_rate, e_rate_ref, max_gap, roh_min_l_initial,
             roh_min_l_final, min_len1, min_len2, cutoff_post, logfile] for iid in iids]
    assert(len(prms[0])==31)   # Sanity Check
                            
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
######################## Functions to estimate Contamination #####################



def hapCon_chrom_BFGS_legacy(iid="", hdf5=None,
                             n_ref=2504, diploid_ref=False, exclude_pops=["AFR"], conPop=["CEU"], 
                             h5_path1000g = None, meta_path_ref = None,folder_out="", c=0.025, roh_jump=300, 
                             e_rate=1e-2, e_rate_ref=1e-3, output=False):
    
    hmm = HMM_Analyze(cython=3, p_model="SardHDF5", e_model="readcount_contam", post_model="Standard",
                      manual_load=True, save=False, save_fp=False, output=output)

    ### Load and prepare the pre-processing Model
    hmm.load_preprocessing_model(conPop)              # Load the preprocessing Model
    hmm.p_obj.set_params(readcounts = True, random_allele=False,
                         folder_out=folder_out, prefix_out_data="", 
                         excluded=exclude_pops, diploid_ref=diploid_ref)
    
    ### Set the paths to ref & target
    hmm.p_obj.set_params(h5_path1000g = h5_path1000g, path_targets = hdf5, 
                         meta_path_ref = meta_path_ref, n_ref=n_ref)
    hmm.load_data(iid=iid, ch='X')  # Load the actual Data
    hmm.load_secondary_objects()
    
    ### Print out the Parameters used in run:
    # print("\nParameters in hapCon_chrom_BFGS:")
    # print("\n".join("{}\t{}".format(k, v) for k, v in parameters.items()))
    # print("\n")

    ### Set the Parameters
    hmm.e_obj.set_params(e_rate = e_rate, e_rate_ref = e_rate_ref)
    hmm.t_obj.set_params(roh_out=0.0, roh_jump=roh_jump)
    con_mle, lower, upper = hmm.optimize_ll_contamination_BFGS(c)
    return con_mle, lower, upper

def hapCon_chrom_BFGS(iid="", mpileup=None, bam=None, bamTable=None, q=30, Q=30,
    n_ref=2504, diploid_ref=False, exclude_pops=["AFR"], conPop=["CEU"], 
    h5_path1000g = None, meta_path_ref = None,
    folder_out="", c=0.025, roh_jump=300, e_rate_ref=1e-3, damage=False,
    logfile=False, output=False, cleanup=False, prefix="hapCon"):
    """Run HapCon to estimate male X chromosome contamination.

    Parameters
    ----------
    iid: str
        IID of the Target Individual, if not provided, will be deduced from the prefix of BAM or mpileup file.
    mpileup: str
        path to mpileup file of chrX.
    bam: str
        path to BAM file.
    bamTable: str
        path to output of BamTable. You need to specify one of mpileup/BAM/bamTable. We recommend using BamTable as your first choice.
    q: int
        minimum mappling quality for reads in BAM file to be considered. Only applicable if used in conjunction with the bam option.
    Q: int
        minimum base quality for reads in BAM file to be considered. Only applicable if used in conjunction with the bam option.
    h5_path1000g: str
        Path of the reference genotypes.
    meta_path_ref: str
        Path of the meta file for the references.
    folder_out: str
        Path of the basis folder for output, if not provided, output will reside in the same directory as BAM or mpileup file.
    output: bool
        Whether to print extensive output. Default False.
    n_ref: int
        Maximum Number of (diploid) reference Individuals to use. Default 2504.
    exclude_pops: list of str
        Which populations to exclude from reference. Default is to exclude African haplotypes.
        Note that the population label supplied in the list must correspond to entries in the "pop" or "super_pop" column in the metadata file. 
        A list of population labels (column "pop") can be seen https://www.coriell.org/1/NHGRI/Collections/1000-Genomes-Collections/1000-Genomes-Project.
        The "super_pop" column has 5 possible values: AFR, EUR, EAS, SAS, AMR.
        The column "pop" is a refinement of "super_pop". For example, to exclude AFR and EUR haplotypes, one can use exclude_pops=[AFR, EUR]. You can also mix "pop" and "super_pop" labels, for example, exluce_pops=[AFR, IBS].
        The same principle applies for the conPop argument as well. 
    conPop: list of str
        use which population in the ref panel as the contaminating pop. If empty list, then use all samples in the ref panel to cauclate allele freq. Default is to CEU allele frequency.
    c: float 
        initial contamination rate to start the BFGS optimization procedure. Default to 0.025.
    roh_jump: float
        Parameter to jump (per Morgan). Default to 300.
    e_rate_ref: float
        Haplotype coping error rate. Default to 0.001.
    logfile: bool
        Whether to produce a logfile. Default False.
    cleanup: bool
        Whether to delete the intermediary .hdf5 file. Default False.
    prefix: str
        The output file will have file name $iid.$prefix.txt
    
    Returns
    ----------
    conMLE: float
        MLE point estimate for contamination.
    lower95: float
        lower bound for the 95% CI of estimated contamination.
    upper95: float
        upper bound for the 95% CI of estimated contamination.
    """    

    if not mpileup and not bam and not bamTable:
        print(f'Must specify one of mpileup/BAM/bamTable.')
        sys.exit()
    elif len(folder_out) == 0:
        if bam:
            folder_out = os.path.dirname(os.path.abspath(bam))
        elif mpileup:
            folder_out = os.path.dirname(os.path.abspath(mpileup))
        else:
            folder_out = os.path.dirname(os.path.abspath(bamTable))

    if len(iid) == 0:
        if bam:
            bamName = os.path.basename(bam)
            iid = bamName[:bamName.find(".bam")]
        elif mpileup:
            mpileupName = os.path.basename(mpileup)
            iid = mpileupName[:mpileupName.find(".mpileup")]
        else:
            bamTableName = os.path.basename(bamTable)
            iid = bamTableName[:bamTableName.find(".BamTable")]
    assert(len(iid) != 0)
    
    # check if the index file for BAM exists
    if bam and not os.path.isfile(bam + ".bai"):
        print(f'No index file found for {bam}. Need to index BAM file using samtools index.')
        sys.exit()

    ### Create Folder if needed, and pipe output if wanted
    prepare_path_general(folder_out, iid, None, "hapCon", logfile) # Set the logfile

    ################## pre-process of mpileup or BAM file ################
    t1 = time.time()
    if bam:
        err, numSitesCovered, path2hdf5 = bam2hdf5(bam, h5_path1000g, ch='X', iid=iid, minMapQual=q, minBaseQual=Q, s=5000000, e=154900000, outPath=folder_out)
        print(f'finished reading bam file, takes {time.time()-t1:.3f}.')
    elif mpileup:
        if not damage:
            err, numSitesCovered, path2hdf5 = mpileup2hdf5(mpileup, h5_path1000g, iid=iid, s=5000000, e=154900000, outPath=folder_out)
        else:
            print(f'Doing damage aware parsing of mpileup file.')
            err, numSitesCovered, path2hdf5 = mpileup2hdf5_damageAware(mpileup, h5_path1000g, iid=iid, s=5000000, e=154900000, outPath=folder_out)
        print(f'finished reading mpileup file, takes {time.time()-t1:.3f}.')
    else:
        err, numSitesCovered, path2hdf5 = bamTable2hdf5(bamTable, h5_path1000g, iid=iid, s=5000000, e=154900000, outPath=folder_out)
        print(f'finished reading BamTable, takes {time.time()-t1:.3f}')

    print(f'number of sites covered by at least one read: {numSitesCovered}')
    print(f'hdf5 file saved to {path2hdf5}')
    if err == 0:
        err = 1e-3 # to avoid numerical errors

    ########################## end of preprocessing ###########################
    # parameters = locals() # Gets dictionary of all local variables at this point
    
    
    hmm = HMM_Analyze(cython=3, p_model="SardHDF5", e_model="readcount_contam", post_model="Standard",
                      manual_load=True, save=False, save_fp=False, output=output)

    ### Load and prepare the pre-processing Model
    hmm.load_preprocessing_model(conPop)              # Load the preprocessing Model
    hmm.p_obj.set_params(readcounts = True, random_allele=False,
                         folder_out=folder_out, prefix_out_data="", 
                         excluded=exclude_pops, diploid_ref=diploid_ref)
    
    ### Set the paths to ref & target
    hmm.p_obj.set_params(h5_path1000g = h5_path1000g, path_targets = path2hdf5, 
                         meta_path_ref = meta_path_ref, n_ref=n_ref)
    hmm.load_data(iid=iid, ch='X')  # Load the actual Data
    hmm.load_secondary_objects()
    
    ### Print out the Parameters used in run:
    # print("\nParameters in hapCon_chrom_BFGS:")
    # print("\n".join("{}\t{}".format(k, v) for k, v in parameters.items()))
    # print("\n")

    ### Set the Parameters
    hmm.e_obj.set_params(e_rate = err/3, e_rate_ref = e_rate_ref)
    hmm.t_obj.set_params(roh_out=0.0, roh_jump=roh_jump)
    #hmm.post_obj.set_params(max_gap=max_gap, cutoff_post=cutoff_post, roh_min_l = roh_min_l)


    con_mle, lower, upper = hmm.optimize_ll_contamination_BFGS(c)
    print(f'estimated contamination rate: {con_mle:.6f}({lower:.6f} - {upper:.6f})')

    #writing output to file
    with open(f'{folder_out}/{iid}.{prefix}.txt', 'w') as out:
        out.write(f'Number of target sites covered by at least one read: {numSitesCovered}\n')
        out.write(f'Method1: Fixing genotyping error rate\n')
        out.write(f'\tEstimated genotyping error via flanking region: {round(err,6)}\n')
        out.write(f'\tMLE for contamination using BFGS: {round(con_mle,6)} ({round(lower,6)} - {round(upper,6)})\n')

    if cleanup:
        os.remove(path2hdf5)
        print(f'delete {path2hdf5}')
    shutil.rmtree(f'{folder_out}/{iid}') # this is a hack

    return con_mle, lower, upper, numSitesCovered



##################################################################################
##################################################################################
################################### END ##########################################


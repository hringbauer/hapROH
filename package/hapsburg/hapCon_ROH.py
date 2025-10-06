# run hapCon_ROH estimate from mpileup results

import numpy as np
import argparse
import time
import sys
import os
from pathlib import Path
import shutil
from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind, hapsb_ind_lowmem
from hapsburg.PackagesSupport.hapsburg_run import hapsb_femaleROHcontam_preload
from hapsburg.PackagesSupport.parallel_runs.helper_functions import multi_run
from hapsburg.PackagesSupport.loadEigenstrat.saveHDF5 import mpileup2hdf5, bamTable2hdf5
from multiprocessing import set_start_method


def main():
    set_start_method("spawn")
    parser = argparse.ArgumentParser(description='Run hapCon_ROH from either mpileup or BamTable output')
    parser.add_argument('--mpileup', action="store", dest="mpath", type=str, required=False,
                        help="Basepath to a list of mpileup file")
    parser.add_argument('--bamTable', action="store", dest="bamTable", type=str, required=False,
                        help="Basepath to a list of BamTable file")          
    parser.add_argument('-i', action="store", dest="iid", type=str, required=True,
                        help="IID of the target individual. Note that the filename of mpileup output or BamTable should follow $iid.chr$ch.mpileup or $iid.chr$ch.BamTable")
    parser.add_argument('-r', action='store', dest='r', type=str, required=True,
                        help='path to reference panel hdf5 file. ')
    parser.add_argument('--meta', action="store", dest="meta", type=str, required=True,
                        help="path to the metadata for the reference panel.")
    parser.add_argument('--minL', action="store", dest="minL", type=float, required=False, default=5.0,
                        help="minimum length of ROH that will be used in contamination estimation (in cM).")
    parser.add_argument('-p', action="store", dest="processes", type=int, required=False, default=1,
                        help="Number of processes to use.")
    parser.add_argument('--niter', action="store", dest="niter", type=int, required=False, default=5,
                        help="Maximum number of iterations.")
    parser.add_argument('--tol', action="store", dest="tol", type=float, required=False, default=0.005,
                        help="Stopping criterion. If the estimated contamination rates between two consecutive iterations differ less than this value, stop iteration.")
    parser.add_argument('--prefix', action="store", dest="prefix", type=str, required=False, default=None,
                        help="prefix of the output. The output will be named as $iid.$prefix.hapCon_ROH.txt")
    parser.add_argument('--lowmem', action="store_true", dest="lowmem", required=False, default=False,
                        help="Use low memory mode.")
    args = parser.parse_args()

    iid = args.iid
    p = args.processes
    mpileup_path= args.mpath
    bamTable_path = args.bamTable
    hapsb_ind_alias = hapsb_ind_lowmem if args.lowmem else hapsb_ind

    if not mpileup_path and not bamTable_path:
        print('Need to at least provide one of pileup files or BamTable files.')
        sys.exit()

    input_path = mpileup_path if mpileup_path else bamTable_path
    abspath = os.path.abspath(input_path)
    basepath = str(Path(abspath).parents[0])
    
    if not os.path.exists(f'{basepath}/hdf5'):
        os.makedirs(f'{basepath}/hdf5')
    print(f'saving hdf5 files in {basepath}/hdf5')

    t1 = time.time()
    if mpileup_path:
        prms = [ [os.path.join(input_path, f'{iid}.chr{ch}.mpileup'),
            args.r + str(ch) + ".hdf5",
            iid, -np.inf, np.inf, f'{basepath}/hdf5', False] for ch in range(1, 23)]
        results = multi_run(mpileup2hdf5, prms, p)
    else:
        prms = [ [os.path.join(input_path, f'{iid}.chr{ch}.BamTable'), 
            args.r + str(ch) + ".hdf5",
            iid, -np.inf, np.inf, f'{basepath}/hdf5', False] for ch in range(1, 23)]
        results = multi_run(bamTable2hdf5, prms, p)
    err = np.mean(np.array([err for err, _, _, _ in results]))
    numSitesCovered = sum([n for _, n, _, _ in results])
    totNumSites = sum([n for _, _, n, _ in results])
    print(f'finished reading input files, takes {round(time.time()-t1, 3)}s')
    print(f'estimated genotyping error: {err}')
    frac = numSitesCovered/totNumSites
    print(f'fraction of sites covered by at least one read: {frac}')
    
    if frac >= 0.7:
        downsample = 1.0
    else:
        downsample = False

    ################################### call ROH ##################################
    ## first, run without contamination
    df = hapsb_ind_alias(iid, chs=range(1,23), 
        path_targets_prefix = f"{basepath}/hdf5",
        h5_path1000g = args.r, meta_path_ref = args.meta,
        folder_out=f"{basepath}/hapRoh_iter/", prefix_out="",
        e_model="readcount_contam", p_model="HDF5", post_model="Standard",
        processes=args.processes, delete=False, output=True, save=True, save_fp=False, 
        n_ref=2504, diploid_ref=True, exclude_pops=[], readcounts=True, random_allele=False, downsample=downsample, 
        c=0.0, roh_min_l_final=0.04, roh_in=1, roh_out=20, roh_jump=300, e_rate=err, e_rate_ref=1e-3, 
        logfile=True, combine=True, file_result="_roh_full.csv")
    df = df[df['lengthM'] >= args.minL/100]
    nBlocks, lengthSum = len(df.index), 100*np.sum(df['lengthM'])
    lengthSum -= 2*nBlocks*0.5
    print(f'number of blocks found with null-model: {nBlocks} with total length {round(lengthSum, 3)}cM')
    # run with 5% contamiantion only if the null-model doesn't yield any ROH or total sum < 10cM
    if nBlocks == 0 or lengthSum < 10:
        hapsb_ind_alias(iid, chs=range(1,23), 
            path_targets_prefix = f"{basepath}/hdf5",
            h5_path1000g = args.r, meta_path_ref = args.meta,
            folder_out=f"{basepath}/hapRoh_iter/", prefix_out="",
            e_model="readcount_contam", p_model="HDF5", post_model="Standard",
            processes=args.processes, delete=False, output=True, save=True, save_fp=False, 
            n_ref=2504, diploid_ref=True, exclude_pops=[], readcounts=True, random_allele=False, downsample=downsample, 
            c=0.05, roh_min_l_final=0.04, roh_in=1, roh_out=20, roh_jump=300, e_rate=err, e_rate_ref=1e-3, 
            logfile=True, combine=True, file_result="_roh_full.csv")
    
    ######################################################################################


    ################## Use the called ROH region to estimate contamination ###############
    contam, se, chunks, sumROH = hapsb_femaleROHcontam_preload(iid, f"{basepath}/hapRoh_iter/{iid}_roh_full.csv",
        args.r, args.meta, hdf5_path=f"{basepath}/hdf5", minLen=args.minL/100,
        e_rate=err, processes=p, prefix=args.prefix, logfile=False, lowmem=args.lowmem)

    # iterate the process if necessary
    if contam >= 0.05:
        print(f'Estimated contamination rate {round(contam, 6)} >0.05 -> Start iteration.')
        diff = np.inf
        niter = 1
        maxIter = args.niter
        contam_prev = contam
        tol = args.tol
        while diff > tol and niter < maxIter:
            hapsb_ind_alias(iid, chs=range(1,23), 
                path_targets_prefix = f"{basepath}/hdf5",
                h5_path1000g = args.r, meta_path_ref = args.meta,
                folder_out=f"{basepath}/hapRoh_iter/", prefix_out="",
                e_model="readcount_contam", p_model="HDF5", post_model="Standard",
                processes=p, delete=True, output=True, save=True, save_fp=False, 
                n_ref=2504, diploid_ref=True, exclude_pops=[], readcounts=True, random_allele=False, downsample=downsample, 
                c=contam_prev, roh_min_l_final=0.05, roh_in=1, roh_out=20, roh_jump=300, e_rate=err, e_rate_ref=1e-3, 
                logfile=True, combine=True, file_result="_roh_full.csv")
            contam, se, chunks, sumROH = hapsb_femaleROHcontam_preload(iid, f"{basepath}/hapRoh_iter/{iid}_roh_full.csv",
                args.r, args.meta, hdf5_path=f"{basepath}/hdf5", minLen=args.minL/100, e_rate=err, 
                processes=p, prefix=args.prefix, logfile=False, lowmem=args.lowmem)
            niter += 1
            diff = abs(contam_prev - contam)
            print(f'iteration {niter} done, prev contam: {round(contam_prev, 6)}, current contam: {round(contam, 6)}')
            contam_prev = contam

        converged = diff <= tol
        if converged:
            print(f'contamination rate converged after {niter} iterations.')
        else:
            print(f'contamination rate did not converge. Try increase the maxIter param.')
        # doing a final round of ROH calling using the new contam estimates
        hapsb_ind_alias(iid, chs=range(1,23), 
                path_targets_prefix = f"{basepath}/hdf5",
                h5_path1000g = args.r, meta_path_ref = args.meta,
                folder_out=f"{basepath}/hapRoh_iter/", prefix_out="",
                e_model="readcount_contam", p_model="HDF5", post_model="Standard",
                processes=p, delete=False, output=True, save=True, save_fp=False, 
                n_ref=2504, diploid_ref=True, exclude_pops=[], readcounts=True, random_allele=False, downsample=downsample,
                c=contam_prev, roh_in=1, roh_out=20, roh_jump=300, e_rate=err, e_rate_ref=1e-3, 
                logfile=True, combine=True, file_result="_roh_full.csv")
    else:
        converged = True


    if args.prefix:
        fileName = f'{iid}.{args.prefix}.hapCon_ROH.txt'
    else:
        fileName = f'{iid}.hapCon_ROH.txt'
    with open(f'{basepath}/{fileName}', 'w') as f:
        f.write(f'Method1: Fixing genotyping error rate at {round(err, 6)}\n')
        f.write(f'\tNumber of ROH blocks found: {len(chunks)}\n')
        f.write(f'\tTotal length of ROH after trimming: {round(100*sumROH,3)}cM\n')
        f.write(f'\tMLE for contamination using BFGS: {round(contam, 6)} ({round(contam-1.96*se, 6)} - {round(contam+1.96*se, 6)})\n')

    
    # clean up
    #shutil.rmtree(f'{basepath}/hdf5') # this is a hack


    #####################################################################################
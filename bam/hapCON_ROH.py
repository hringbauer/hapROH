# run hapCON on a known list of ROH blocks.
# these ROH blocks will be treated as ground-truth, so no iterative procedure will be performed.

import argparse
import sys
import time
import os
import numpy as np
from toHDF5 import mpileup2hdf5



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run hapROH(on autosomes) from a list of known ROH blocks')
    parser.add_argument('--hdf5', action="store", dest="hdf5", type=str, required=False,
                        help="Basepath to a list of hdf5 files")
    parser.add_argument('--mpileup', action="store", dest="mpileup", type=str, required=False,
                        help="Basepath to a list of mpileup file")
    parser.add_argument('-r', action="store", dest="roh", type=str, required=True,
                        help='file containing a list of known ROH blocks.')
    parser.add_argument('-e', action="store", dest="err", type=float, required=False,
                        help='genotyping error rate')
    parser.add_argument('-i', action="store", dest="iid", type=str, required=True,
                        help="IID of the target individual.")
    parser.add_argument('-p', action="store", dest="processes", type=int, required=False, default=1,
                        help="Number of processes to use.")
    args = parser.parse_args()

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")
    from hapsburg.PackagesSupport.hapsburg_run import hapsb_femaleROHcontam_preload
    from hapsburg.PackagesSupport.parallel_runs.helper_functions import multi_run
    from multiprocessing import set_start_method
    set_start_method("spawn")

    iid = args.iid
    p = args.processes
    hdf5_path= args.hdf5
    roh_path = args.roh
    err = args.err

    if hdf5_path == None:
        mpileup_path= args.mpileup
        if mpileup_path == None:
            print(f'need to at least provide a path to hdf5 files or a path fo mpileup files...')
            sys.exit()
        basepath = ""
        if not mpileup_path.endswith('/'):
            basepath = mpileup_path[:mpileup_path.rindex("/")]
            mpileup_path += '/'
        else:
            basepath = mpileup_path[:mpileup_path.rindex("/", end=len(mpileup_path)-1)]

        hdf5_path = f'{basepath}/hdf5'
        if not os.path.exists(hdf5_path):
            os.makedirs(hdf5_path)
            print(f'saving hdf5 files in {hdf5_path}')

        t1 = time.time()
        prms = [ [mpileup_path + f'{iid}.chr{ch}.mpileup', 
            f'/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr{ch}.hdf5',
            iid, -np.inf, np.inf, f'{basepath}/hdf5', False] for ch in range(1, 23)]
        results = multi_run(mpileup2hdf5, prms, p)
        err = np.mean(np.array([err for err, _, _ in results]))
        numSitesCovered = sum([n for _, n, _ in results])
        print(f'finished reading mpileup files, takes {round(time.time()-t1, 3)}s')
        print(f'estimated genotyping error: {err}')
    else:
        basepath = ""
        if not hdf5_path.endswith('/'):
            basepath = hdf5_path[:hdf5_path.rindex("/")]
            hdf5_path += '/'
        else:
            basepath = hdf5_path[:hdf5_path.rindex("/", end=len(hdf5_path)-1)]

    contam, se = hapsb_femaleROHcontam_preload(iid, roh_path, hdf5_path,
        "/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr", 
        "/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv",
        f'{basepath}', init_c=0.025, conPop=["CEU"],
        roh_jump=300, e_rate=err, e_rate_ref=1e-3,
        processes=p, n_ref=2504, diploid_ref=True, 
        exclude_pops=[], p_model="SardHDF5", logfile=True)
    
    # with open(f'{basepath}/{iid}.MyRoh.1240k.results', 'w') as f:
    #     f.write(f'Method1: Fixing genotyping error rate at {err}\n')
    #     f.write(f'\tROH blocks obtained from: {roh_path}\n')
    #     f.write(f'\tMLE for contamination using BFGS: {round(contam, 6)} ({round(contam-1.96*se, 6)} - {round(contam+1.96*se, 6)})\n')
    
    


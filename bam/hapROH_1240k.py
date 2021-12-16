import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os as os
import sys as sys
import multiprocessing as mp
import argparse
import time
from toHDF5 import mpileup2hdf5



## Call ROH
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run hapROH(on autosomes) from mpileup output')
    parser.add_argument('-b', action="store", dest="mpath", type=str, required=False,
                        help="Basepath to a list of mpileup file")
    parser.add_argument('-i', action="store", dest="iid", type=str, required=True,
                        help="IID of the target individual.")
    parser.add_argument('-p', action="store", dest="processes", type=int, required=False, default=1,
                        help="Number of processes to use.")
    
    args = parser.parse_args()

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")
    from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind
    from hapsburg.PackagesSupport.parallel_runs.helper_functions import multi_run
    from multiprocessing import set_start_method
    set_start_method("spawn")

    iid = args.iid
    p = args.processes
    mpileup_path= args.mpath

    basepath = ""
    if not mpileup_path.endswith('/'):
        basepath = mpileup_path[:mpileup_path.rindex("/")]
        mpileup_path += '/'
    else:
        basepath = mpileup_path[:mpileup_path.rindex("/", end=len(mpileup_path)-1)]
    
    if not os.path.exists(f'{basepath}/hdf5'):
        os.makedirs(f'{basepath}/hdf5')
    print(f'saving hdf5 files in {basepath}/hdf5')

    t1 = time.time()
    prms = [ [mpileup_path + f'{iid}.chr{ch}.mpileup', 
            f'/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr{ch}.hdf5',
            iid, -np.inf, np.inf, f'{basepath}/hdf5', False] for ch in range(1, 23)]
    results = multi_run(mpileup2hdf5, prms, p)
    err = np.mean(np.array([err for err, _, _ in results]))
    numSitesCovered = sum([n for _, n, _ in results])
    print(f'finished reading mpileup files, takes {round(time.time()-t1, 3)}s')
    print(f'estimated genotyping error: {err}')

    ################################### call ROH ##################################
    hapsb_ind(iid, chs=range(1,23), 
        path_targets_prefix = f"{basepath}/hdf5",
        h5_path1000g = "/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr", 
        meta_path_ref = "/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv",
        folder_out=f"{basepath}/hapRoh/", prefix_out="",
        e_model="readcount", p_model="SardHDF5", post_model="Standard",
        processes=args.processes, delete=False, output=True, save=True, save_fp=False, 
        n_ref=2504, diploid_ref=True, exclude_pops=[], readcounts=True, random_allele=False,
        e_rate=err, e_rate_ref=1e-3, logfile=True, combine=True, file_result="_roh_full.csv")
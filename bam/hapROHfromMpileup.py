# run hapROH from mpileup results

import numpy as np
import argparse
import time
import sys
from toHDF5 import mpileup2hdf5


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run hapROH(on autosomes) from mpileup output')
    parser.add_argument('-b', action="store", dest="basepath", type=str, required=False,
                        help="Basepath to a list of mpileup file")
    parser.add_argument('-i', action="store", dest="iid", type=str, required=True,
                        help="IID of the target individual.")
    parser.add_argument('-p', action="store", dest="processes", type=int, required=False, default=1,
                        help="Number of processes to use.")
    args = parser.parse_args()

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")
    from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind
    from hapsburg.PackagesSupport.parallel_runs.helper_functions import multi_run

    iid = args.iid
    p = args.processes
    basepath= args.basepath
    if not basepath.endswith('/'):
        basepath += '/'
    
    t1 = time.time()
    prms = [ [basepath + f'{iid}.chr{ch}.mpileup', 
            f'/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr{ch}.hdf5',
            iid, -np.inf, np.inf, './hdf5', False] for ch in range(1, 23)]
    results = multi_run(mpileup2hdf5, prms, p)
    err = np.mean(np.array([err for err, _, _ in results]))
    print(f'finished reading mpileup file, takes {round(time.time()-t1, 3)}s')
    print(f'estimated genotyping error: {err}')

    # for ch in range(1, 23):
    #     t1 = time.time()
    #     path2mpileup = args.basepath
    #     if not path2mpileup.endswith('/'):
    #         path2mpileup += "/"
    #     path2mpileup += f'{iid}.chr{ch}.mpileup'
    #     path2ref = f'/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr{ch}.hdf5'
    #     err, numSitesCovered, path2hdf5 = mpileup2hdf5(path2mpileup, path2ref, iid=iid, outPath='./hdf5', output=False)
    #     print(f'finished reading mpileup file for chr{ch}, takes {round(time.time()-t1, 3)}.')
    #     print(f'number of sites covered by at least one read: {numSitesCovered}')
    #     print(f'hdf5 file saved to {path2hdf5}')

    # hapsb_ind(iid, chs=range(1,23), 
    #     path_targets_prefix = f"/mnt/archgen/users/yilei/Data/AGDP/contamX/{iid}/hdf5",
    #     h5_path1000g = "/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr", 
    #     meta_path_ref = "/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv",
    #     folder_out=f"/mnt/archgen/users/yilei/Data/AGDP/contamX/{iid}/hapRoh/", prefix_out="",
    #     e_model="readcount", p_model="SardHDF5", post_model="Standard",
    #     processes=args.threads, delete=False, output=True, save=True, save_fp=False, 
    #     n_ref=2504, diploid_ref=True, exclude_pops=[], readcounts=True, random_allele=False,
    #     roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=1e-3, 
    #     cutoff_post = 0.999, max_gap=0, roh_min_l = 0.04, logfile=True, combine=True, 
    #     file_result="_roh_full.csv")
    

# run hapCon on single sample across all autosomes (for male ones, chrX too but this will be added later)

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os as os
import sys as sys
import multiprocessing as mp


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='run hapCon on single sample across all autosomes. For male ones, chrX too but this will be added later')
    parser.add_argument('-i', action="store", dest="iid", type=str, required=True, 
                        help='individual id as found in the target file')
    args = parser.parse_args()

    print(f"CPU Count: {mp.cpu_count()}")
    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")
    #path = "/mnt/archgen/users/yilei"   # The Path on Yilei's remote space
    #os.chdir(path)
    #print(f"Set path to: {os.getcwd()}") # Show the current working directory.

    from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom  # Need this import
    from hapsburg.PackagesSupport.hapsburg_run import hapCon_ind

    outPath = "/mnt/archgen/users/yilei/tools/hapROH/test_output/contaminate/testScript/"

    loglls = []
    cons = np.arange(0, 0.2, 0.005)
    for con in cons:
        ll = hapCon_ind(args.iid, chs=range(1,23), 
            path_targets='/mnt/archgen/users/yilei/Data/SA_1240KHDF5/marcus2020.h5',
            h5_path1000g='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr',
            meta_path_ref='/mnt/archgen/users/yilei/Data/1000G/meta_df_all.csv', 
            folder_out=outPath, prefix_out="",
            e_model="readcount_contam", p_model="SardHDF5", post_model="Standard",
            processes=1, delete=False, output=False, save=True, save_fp=False, 
            n_ref=2504, diploid_ref=True, exclude_pops=[], readcounts=True, random_allele=False,
            c=con, roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.00, 
            cutoff_post = 0.999, max_gap=0, roh_min_l = 0.01, logfile=False)
        loglls.append(ll)
    
    # plotting loglikelihood versus con rate
    plt.plot(cons, loglls)
    plt.xlabel('contamination rate')
    plt.ylabel('loglikelihoods')
    plt.title(args.iid)
    plt.savefig(f'{outPath}/{args.iid}.png', dpi=300)









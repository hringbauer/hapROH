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
    parser.add_argument('-t', action="store", dest="target", type=str, required=False,
                        default="/mnt/archgen/users/yilei/Data/SA_1240KHDF5/marcus2020.h5",
                        help="path to HDF5 file for the target individual")
    parser.add_argument('-o', action="store", dest="o", type=str, required=False,
                        default="/mnt/archgen/users/yilei/tools/hapROH/test_output/contaminate/testScript/",
                        help="output path")
    parser.add_argument('--chr', action="store", dest="chr", type=str, required=False,
                        help="chromosomes to run hapCon on. In the format 1-22,X for example.")
    parser.add_argument('--er', action="store", dest="e_rate_ref", type=float, required=False, default=1e-3,
                        help="Reference panel error rate. Used to emulate mutation since common ancestry with the copied haplotype.")
    args = parser.parse_args()

    print(f"CPU Count: {mp.cpu_count()}")
    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")

    from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom  # Need this import
    from hapsburg.PackagesSupport.hapsburg_run import hapCon_ind

    # parsing chromosome string
    chs = []
    if len(args.chr) == 0:
        chs.extend(range(1, 23))
    else:
        for item in args.chr.split(','):
            i = item.find('-')
            if i != -1:
                start, end = int(item[:i]), int(item[i+1:])
                chs.extend(range(start, 1+end))
            else:
                if item in ['X', 'x']:
                    chs.append(23)
                else:
                    chs.append(int(item))
    print(f'running hapCon on chromosome {chs}')

    cons = np.arange(0, 0.2, 0.005)
    lls, con_mle, lower, upper = hapCon_chrom(args.iid, ch='X', 
            path_targets=args.target,
            h5_path1000g='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr',
            meta_path_ref='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv', 
            folder_out=args.o, prefix_out="fast_",
            e_model="readcount_contam", p_model="SardHDF5", post_model="Standard", save=True, save_fp=False, 
            n_ref=2504, diploid_ref=True, exclude_pops=["TSI"], conPop=[], readcounts=True, random_allele=False,
            c=cons, roh_in=1, roh_out=0, roh_jump=300, e_rate=0.001, e_rate_ref=args.e_rate_ref, 
            cutoff_post = 0.999, max_gap=0, roh_min_l = 0.01, logfile=False)
    
    # plotting loglikelihood versus con rate
    plt.plot(cons, lls)
    plt.xlabel('contamination rate')
    plt.ylabel('loglikelihoods')
    plt.axvline(con_mle, linestyle='dashed')
    plt.axvline(lower, linestyle='dashed', color='red')
    plt.axvline(upper, linestyle='dashed', color='red')
    plt.title(f'{args.iid}\t{con_mle}')
    plt.savefig(f'{args.o}/{args.iid}.png', dpi=300)









import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os as os
import sys as sys
import multiprocessing as mp
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run hapCon on Mosaic X chromosome with readcount data.')
    parser.add_argument('--cov', action="store", dest="cov", type=float, required=True,
                        help="Genomic coverage")
    parser.add_argument('--con', action="store", dest="con", type=float, required=True,
                        help="Contamination Rate")
    parser.add_argument('--err', action="store", dest="err", type=float, required=True,
                        help="genotyping error.")
    parser.add_argument('--eref', action="store", dest="eref", type=float, required=True,
                        help="error rate when copied from the reference panel.")
    parser.add_argument('--conpop1', action="store", dest="conpop1", type=str, required=True,
                        help="which source of population is used as contamination in simulation")
    parser.add_argument('--conpop2', action="store", dest="conpop2", type=str, required=True,
                        help="which source of population is used as contamination in inference")
    args = parser.parse_args()

    path = "/mnt/archgen/users/yilei/tools/hapROH"   # The Path on Yilei's remote space
    os.chdir(path)  # Set the right Path (in line with Atom default)

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")  # hack to get local package first in path [FROM HARALD - DELETE!!!]
    from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom  # Need this import
    from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom_BFGS  # Need this import


    base_path="./simulated/1000G_Mosaic/TSI/maleXMisAnc/" 
    path1000G="/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr"
    ch='X'

    # parameters for readcount data  
    cov = args.cov
    con = args.con
    err_rate = args.err
    e_rate_ref = args.eref

    base_path += f'{args.conpop1}/'

    if con == 0:
        base_path += "con0/"
    elif con == 0.05:
        base_path += "con5/"
    elif con == 0.1:
        base_path += "con10/"
    
    prefix = ""
    if cov == 0.05:
        prefix = "chrX_cov1over20"
    elif cov == 0.1:
        prefix = "chrX_cov1over10"
    elif cov == 0.5:
        prefix = "chrX_cov1over2"
    elif cov == 1.0:
        prefix = "chrX_cov1"
    elif cov == 2.0:
        prefix = "chrX_cov2"
    elif cov == 5.0:
        prefix = "chrX_cov5"

    outFolder = base_path + prefix

    results = np.zeros((100, 3))
    for i in range(100):
        iid = "iid" + str(i)
        conMLE, lower95, upper95 = hapCon_chrom_BFGS(iid, ch='X', save=False, save_fp=False, 
            n_ref=2504, diploid_ref=True, exclude_pops=["PEL"], conPop=[args.conpop2], 
            e_model="readcount_contam", p_model="SardHDF5", readcounts=True, random_allele=False,
            post_model="Standard", 
            path_targets=f"{outFolder}/data.h5",
            h5_path1000g='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr',
            meta_path_ref='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv', 
            folder_out=outFolder, prefix_out="",
            c=0.025, roh_in=1, roh_out=0, roh_jump=300, e_rate=err_rate, e_rate_ref=e_rate_ref,
            max_gap=0, cutoff_post = 0.999, roh_min_l = 0.01, logfile=False)

        results[i, :] = (conMLE, lower95, upper95)
    
    # write output to a file
    with open(f'{outFolder}/batchresults_bfgs_{args.conpop1}_{args.conpop2}.txt', 'w') as out:
        out.write(f'###contamination={con}, coverage={cov}, genotyping error={err_rate}, ref err={e_rate_ref}\n')
        out.write(f'###sampleID\tconMLE\tlower95CI\tupper95CI\n')
        for i in range(100):
            iid = "iid" + str(i)
            conMLE, lower95, upper95 = results[i]
            out.write(f'{iid}\t{conMLE}\t{lower95}\t{upper95}\n')
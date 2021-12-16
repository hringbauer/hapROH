import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os as os
import sys as sys
import multiprocessing as mp
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run hapCon on Mosaic X chromosome with readcount data. Estimate contammination and genotyping error rate jointly')
    parser.add_argument('--cov', action="store", dest="cov", type=float, required=True,
                        help="Genomic coverage")
    parser.add_argument('--con', action="store", dest="con", type=float, required=True,
                        help="Contamination Rate")
    parser.add_argument('--err', action="store", dest="err", type=float, required=False, default=1e-2,
                        help="genotyping error.")
    parser.add_argument('--eref', action="store", dest="eref", type=float, required=True,
                        help="error rate when copied from the reference panel.")
    parser.add_argument('-b', action="store", dest="base", type=str, required=True, help="base path")
    args = parser.parse_args()

    path = "/mnt/archgen/users/yilei/tools/hapROH"   # The Path on Yilei's remote space
    os.chdir(path)  # Set the right Path (in line with Atom default)

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")  # hack to get local package first in path [FROM HARALD - DELETE!!!]
    from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom_2d  # Need this import
    
    #base_path="./simulated/1000G_Mosaic/TSI/maleX7/" 
    base_path = args.base
    path1000G="/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/chX1240/chr"
    ch='X'

    # parameters for readcount data  
    cov = args.cov
    con = args.con
    err_rate = args.err
    e_rate_ref = args.eref

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

    results = np.zeros((100, 6))
    for i in range(100):
        iid = "iid" + str(i)
        mle, se = hapCon_chrom_2d(iid, ch='X', 
            path_targets=f"{outFolder}/data.h5",
            h5_path1000g='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr',
            meta_path_ref='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv', 
            folder_out=outFolder, prefix_out="",
            e_model="readcount_contam", p_model="SardHDF5", post_model="Standard", save=False, save_fp=False, 
            n_ref=2504, diploid_ref=True, exclude_pops=["TSI"], conPop=[], readcounts=True, random_allele=False,
            c=0.025, roh_in=1, roh_out=0, roh_jump=300, e_rate=0.01, e_rate_ref=1e-3, 
            cutoff_post = 0.999, max_gap=0, roh_min_l = 0.01, logfile=False)
        results[i, :] = (mle[0], mle[0]-1.96*se[0], mle[0]+1.96*se[0], mle[1], mle[1]-1.96*se[1], mle[1]+1.96*se[1])
    
    # write output to a file
    with open(f'{outFolder}/batchresults_2d.txt', 'w') as out:
        out.write(f'###contamination={con}, coverage={cov}, ref err={e_rate_ref}\n')
        out.write(f'###sampleID\tmle_contamination\tcon_low95\tcon_high95\tmle_err\terr_low95\terr_high95\n')
        for i in range(100):
            iid = "iid" + str(i)
            a, b, c, d, e, f = results[i]
            out.write(f'{iid}\t{a}\t{b}\t{c}\t{d}\t{e}\t{f}\n')
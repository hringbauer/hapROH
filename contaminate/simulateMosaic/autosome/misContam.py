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
    parser.add_argument('--conI', action="store", dest="conI", type=float, required=True,
                        help="Contamination Rate used for inference")
    parser.add_argument('--conS', action="store", dest="conS", type=float, required=True,
                        help="Contamination Rate used in simulation")
    parser.add_argument('--nblock', action="store", dest="nblock", type=int, required=True,
                        help="number of ROH blocks")
    parser.add_argument('-c', action="store", dest="chr", type=int, required=False, default=1,
                        help="Which autosome to decode. Default is chr1.")
    parser.add_argument('--err', action="store", dest="err", type=float, required=False, default=1e-2,
                        help="genotyping error.")
    parser.add_argument('--eref', action="store", dest="eref", type=float, required=False, default=1e-3,
                        help="error rate when copied from the reference panel.")
    args = parser.parse_args()

    path = "/mnt/archgen/users/yilei/tools/hapROH"   # The Path on Yilei's remote space
    os.chdir(path)  # Set the right Path (in line with Atom default)

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")  # hack to get local package first in path [FROM HARALD - DELETE!!!]
    from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind # Need this import


    base_path="./simulated/1000G_Mosaic/CHB/Autosome/" 
    path1000G="/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr"
    ch=args.chr

    # parameters for readcount data  
    cov = args.cov
    conS = args.conS
    conI = args.conI
    err_rate = args.err
    e_rate_ref = args.eref

    nblocks = args.nblock

    if conS == 0:
        base_path += "con0/"
    elif conS == 0.05:
        base_path += "con5/"
    elif conS == 0.1:
        base_path += "con10/"
    elif conS == 0.15:
        base_path += "con15/"
    elif conS == 0.2:
        base_path += "con20/"
    elif conS == 0.25:
        base_path += "con25/"
    
    base_path += f'{nblocks}blocks/'
    
    prefix = ""
    if cov == 0.05:
        prefix = "cov1over20"
    elif cov == 0.1:
        prefix = "cov1over10"
    elif cov == 0.5:
        prefix = "cov1over2"
    elif cov == 1.0:
        prefix = "cov1"
    elif cov == 2.0:
        prefix = "cov2"
    elif cov == 5.0:
        prefix = "cov5"

    targetFolder = base_path + prefix

    # prepare output folder path
    outFolder = "./simulated/1000G_Mosaic/CHB/misContam/"
    if conS == 0:
        outFolder += "conS0/"
    elif conS == 0.05:
        outFolder += "conS5/"
    elif conS == 0.1:
        outFolder += "conS10/"
    elif conS == 0.15:
        outFolder += "conS15/"
    elif conS == 0.2:
        outFolder += "conS20/"
    elif conS == 0.25:
        outFolder += "conS25/"

    if conI == 0:
        outFolder += "conI0/"
    elif conI == 0.02:
        outFolder += "conI2/"
    elif conI == 0.05:
        outFolder += "conI5/"
    elif conI == 0.075:
        outFolder += "conI75/"
    elif conI == 0.1:
        outFolder += "conI10/"
    elif conI == 0.15:
        outFolder += "conI15/"
    elif conI == 0.2:
        outFolder += "conI20/"
    elif conI == 0.25:
        outFolder += "conI25/"

    if cov == 0.5:
        outFolder += "cov1over2"
    elif cov == 1.0:
        outFolder += "cov1"
    elif cov == 2.0:
        outFolder += "cov2"
    elif cov == 5.0:
        outFolder += "cov5"


    for i in range(100):
        iid = "iid" + str(i)
        hapsb_ind(iid, chs=range(1,2), 
        path_targets = f"{targetFolder}/data.h5",
        h5_path1000g = "/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr", 
        meta_path_ref = "/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv",
        folder_out=f"{outFolder}/hapRoh/", prefix_out="",
        e_model="readcount_contam", p_model="SardHDF5", post_model="Standard",
        processes=1, delete=True, output=True, save=True, save_fp=False, 
        c=conI, conPop=["CEU"],
        n_ref=2504, diploid_ref=True, exclude_pops=["CHB"], readcounts=True, random_allele=False,
        roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=1e-3, 
        cutoff_post = 0.999, max_gap=0.005, roh_min_l = 0.04, logfile=False, combine=True, 
        file_result="_roh_full.csv")
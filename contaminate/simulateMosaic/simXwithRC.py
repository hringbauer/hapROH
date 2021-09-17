import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os as os
import sys as sys
import multiprocessing as mp
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate Mosaic X chromosome with readcount data.')
    parser.add_argument('-n', action="store", dest="n", type=int, required=True, 
                        help='Number of haplotypes to simulate.')
    parser.add_argument('--cov', action="store", dest="cov", type=float, required=True,
                        help="Genomic coverage")
    parser.add_argument('--con', action="store", dest="con", type=float, required=True,
                        help="Contamination Rate")
    parser.add_argument('--err', action="store", dest="err", type=float, required=True,
                        help="genotyping error.")
    parser.add_argument('--eref', action="store", dest="eref", type=float, required=True,
                        help="error rate when copied from the reference panel.")
    parser.add_argument('--prefix', action="store", dest="prefix", type=str, required=False, default="",
                        help="Prefix for the output file.")
    parser.add_argument('-b', action="store", dest="b", type=str, required=False,
                        default="/mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleX/",
                        help="output base path")
    args = parser.parse_args()

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/create1000G_Mosaic")
    from createMosaicsMulti import create_individual_mosaic

    path = "/mnt/archgen/users/yilei/tools/hapROH"   # The Path on Yilei's remote space
    os.chdir(path)  # Set the right Path (in line with Atom default)
    
    # actual simulation
    base_path=args.b
    path1000G="/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/chX1240/chr"
    pop_list=["TSI"]
    conPop=[]
    n = args.n  # Number of Individuals
    ch='X'
    chunk_length=0.0025
    l = 5  # length of simulated ROH, in centiMorgen
    n_blocks=5 # How many blocks will be copied in

    # parameters for faking readcount data  
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
    
    if cov == 0.05:
        base_path += "chrX_cov1over20/"
    elif cov == 0.1:
        base_path += "chrX_cov1over10/"
    elif cov == 0.5:
        base_path += "chrX_cov1over2/"
    elif cov == 1.0:
        base_path += "chrX_cov1/"
    elif cov == 2.0:
        base_path += "chrX_cov2/"
    elif cov == 5.0:
        base_path += "chrX_cov5/"

    print(f'simulating {n} maleX chromosomes')
    print(f'output basepath: {base_path}')
    create_individual_mosaic(base_path, path1000G, pop_list, n, ch, chunk_length, l, n_blocks, cov, con, err_rate, e_rate_ref, conPop, args.prefix)

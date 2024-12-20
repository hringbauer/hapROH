# simulate one autosome (with contamination or not)
import numpy as np
import pandas as pd
import os as os
import sys as sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate Mosaic autosome with readcount data.')
    parser.add_argument('-n', action="store", dest="n", type=int, required=True, 
                        help='Number of haplotypes to simulate.')
    parser.add_argument('--cov', action="store", dest="cov", type=float, required=True,
                        help="Genomic coverage")
    parser.add_argument('--con', action="store", dest="con", type=float, required=True,
                        help="Contamination Rate")
    parser.add_argument('--err', action="store", dest="err", type=float, required=False, default=1e-2,
                        help="genotyping error.")
    parser.add_argument('--eref', action="store", dest="eref", type=float, required=False, default=1e-3,
                        help="error rate when copied from the reference panel.")
    parser.add_argument('--nblock', action="store", dest="nblock", type=int, required=True,
                        help="number of ROH blocks to simulate.")
    parser.add_argument('-l', action="store", dest="length", type=float, required=False, default=8.0,
                        help="length of each ROH block, default to 8.0cM.")
    parser.add_argument('-c', action="store", dest="chr", type=int, required=False, default=1,
                        help="Which autosome to simulate. Default is chr1.")
    parser.add_argument('--downsample', action="store_true", help="whether to downsample numer of reads at each site to 1.")
    parser.add_argument('--hetero', action="store_true", help="whether to simulate uneven coverage of 1240k targets.")

    args = parser.parse_args()

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/create1000G_Mosaic")
    from createMosaicsMulti import create_individual_mosaic

    path = "/mnt/archgen/users/yilei/tools/hapROH"   # The Path on Yilei's remote space
    os.chdir(path)  # Set the right Path (in line with Atom default)
    
    # actual simulation
    base_path='/mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/CHB/Autosome_wgs/'
    path1000G="/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/maf5_auto/maf5_chr"
    pop_list=["CHB"]
    conPop=["CEU"]
    n = args.n  # Number of Individuals
    ch=args.chr
    chunk_length=0.0025
    l =  args.length # length of simulated ROH, in centiMorgen
    n_blocks=args.nblock # How many blocks will be copied in

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
    elif con == 0.15:
        base_path += "con15/"
    elif con == 0.2:
        base_path += "con20/"
    elif con == 0.25:
        base_path += "con25/"
    
    base_path += f'{n_blocks}blocks/'
    if cov == 0.05:
        base_path += "cov1over20/"
    elif cov == 0.1:
        base_path += "cov1over10/"
    elif cov == 0.5:
        base_path += "cov1over2/"
    elif cov == 1.0:
        base_path += "cov1/"
    elif cov == 2.0:
        base_path += "cov2/"
    elif cov == 5.0:
        base_path += "cov5/"

    print(f'simulating {n} chr{ch}. Each with {n_blocks} length {l}cM ROH blocks.')
    print(f'output basepath: {base_path}')
    create_individual_mosaic(base_path=base_path, path1000G=path1000G, pop_list=pop_list, 
        n=n, ch=ch, l=l, n_blocks=n_blocks, cov=cov, con=con, err_rate=err_rate, e_rate_ref=e_rate_ref, 
        conPop=conPop, downsample=args.downsample, heterogeneous=args.hetero)


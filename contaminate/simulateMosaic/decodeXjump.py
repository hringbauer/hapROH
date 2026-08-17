import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os as os
import sys as sys
import multiprocessing as mp
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run hapCon on Mosaic X chromosome with readcount data.')
    parser.add_argument('-j', action="store", dest="jump", type=int, required=True,
                        help="jump rate index (see simulation readme).")
    parser.add_argument('--cov', action="store", dest="cov", type=float, required=True,
                        help="Genomic coverage")
    args = parser.parse_args()

    path = "/mnt/archgen/users/yilei/tools/hapROH"   # The Path on Yilei's remote space
    os.chdir(path)  # Set the right Path (in line with Atom default)

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")  # hack to get local package first in path [FROM HARALD - DELETE!!!]
    from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom_BFGS_legacy  # Need this import
    
    base_path="./simulated/1000G_Mosaic/TSI/maleXjump/" 
    path1000G="/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/chX1240/chr"
    ch='X'

    # parameters for readcount data  
    cov = args.cov
    err_rate = 1e-2
    e_rate_ref = 1e-3
    jump = args.jump
    
    prefix = ""
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

    outFolder = base_path + "jump" + str(jump)
    os.system(f'rm -r {outFolder}/iid*')


    results = np.zeros((100, 3))
    for i in range(100):
        iid = "iid" + str(i)
        conMLE, lower95, upper95 = hapCon_chrom_BFGS_legacy(iid,
            hdf5=f"{outFolder}/data.h5",
            h5_path1000g='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/chr',
            meta_path_ref='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv', 
            folder_out=outFolder, exclude_pops=["TSI"], conPop=[], \
            roh_jump=300, e_rate=err_rate, e_rate_ref=e_rate_ref)
        results[i, :] = (conMLE, lower95, upper95)
    
    # write output to a file
    with open(f'{outFolder}/batchresults.txt', 'w') as out:
        out.write(f'###contamination=0.1, coverage={cov}, copying jump rate index={jump}, genotyping error used in inference=1e-2, ref err=1e-3\n')
        out.write(f'###sampleID\tconMLE\tlower95CI\tupper95CI\n')
        for i in range(100):
            iid = "iid" + str(i)
            conMLE, lower95, upper95 = results[i]
            out.write(f'{iid}\t{conMLE}\t{lower95}\t{upper95}\n')
            
    os.system(f'rm -r {outFolder}/iid*')

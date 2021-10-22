import argparse
import time
import os
import sys
from toHDF5 import mpileup2hdf5, bam2hdf5


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert bam file to hdf5 format that stores readcount info at target sites.')
    parser.add_argument('-m', action="store", dest="mpileup", type=str, required=False,
                        help="path to samtools mpileup file")
    parser.add_argument('-b', action="store", dest="bam", type=str, required=False,
                        help="path to bam file")
    parser.add_argument('-r', action="store", dest="ref", type=str, required=True,
                        help="path to reference panel")
    parser.add_argument('-o', action="store", dest="out", type=str, required=False, default="",
                        help="path to the hdf5 file output")
    parser.add_argument('-i', action="store", dest="iid", type=str, required=False, default="",
                        help="IID of the target individual. If unspecified, will use the prefix of the bam file.")
    parser.add_argument('-t', action="store", dest="trim", required=False, default=0, 
                        help="trim certain number of bases from both ends.")
    args = parser.parse_args()

    iid = args.iid
    if len(iid) == 0:
        if args.bam != None:
            bamName = os.path.basename(args.bam)
            iid = bamName[:bamName.find(".bam")]
        elif args.mpileup != None:
            mpileupName = os.path.basename(args.mpileup)
            iid = mpileupName[:mpileupName.find(".mpileup")]
    assert(len(iid) != 0)

    t1 = time.time()
    if args.mpileup != None:
        err, numSitesCovered, path2hdf5 = mpileup2hdf5(args.mpileup, args.ref, iid=iid, s=5000000, e=154900000, outPath=args.out)
        print(f'finished reading mpileup file, takes {round(time.time()-t1, 3)}.')
    if args.bam != None:
        err, numSitesCovered, path2hdf5 = bam2hdf5(args.bam, args.ref, iid=iid, s=5000000, e=154900000, outPath=args.out, trim=args.trim)
        print(f'finished reading bam file, takes {round(time.time()-t1, 3)}.')
        
    print(f'number of sites covered by at least one read: {numSitesCovered}')
    print(f'hdf5 file saved to {path2hdf5}')

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")
    from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom_BFGS

    #err = err/3.0

    mle_bfgs, low95_bfgs, up95_bfgs = hapCon_chrom_BFGS(iid, 'X', save=False, save_fp=False, n_ref=2504, diploid_ref=False, 
        exclude_pops=[], conPop=["CEU"], e_model="readcount_contam", p_model="SardHDF5", 
        readcounts=True, random_allele=False, post_model="Standard", path_targets=path2hdf5, 
        folder_out='/mnt/archgen/users/yilei/Data/iberian_BAM/hapCon/', h5_path1000g=args.ref,
        meta_path_ref='/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv', 
        prefix_out=iid, c=0.025, roh_in=1, roh_out=0, roh_jump=300, e_rate=err, e_rate_ref=1e-3,
        max_gap=0, cutoff_post = 0.999, roh_min_l = 0.01, logfile=False)

    with open(f'{iid}.hapcon.CEU.txt', 'w') as out:
        out.write(f'Number of target sites covered by at least one read: {numSitesCovered}\n')
        out.write(f'Method1: Fixing genotyping error rate\n')
        out.write(f'\tEstimated genotyping error via flanking region: {round(err, 6)}\n')
        out.write(f'\tMLE for contamination using BFGS: {round(mle_bfgs, 6)} ({round(low95_bfgs, 6)} - {round(up95_bfgs, 6)})\n')
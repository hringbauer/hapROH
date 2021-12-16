import argparse
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert bam file to hdf5 format that stores readcount info at target sites.')
    parser.add_argument('-m', action="store", dest="mpileup", type=str, required=False,
                        help="path to samtools mpileup file")
    parser.add_argument('-b', action="store", dest="bam", type=str, required=False,
                        help="path to bam file")
    parser.add_argument('-r', action="store", dest="ref", type=str, required=True,
                        help="path to reference panel")
    parser.add_argument('--meta', action="store", dest="meta", type=str, required=False, 
                        default="/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv",
                        help="path to metadata of the reference panel")
    parser.add_argument('--con', action="store", dest="conpop", type=str, required=False,
                        help='source of contamination. Default is CEU.')
    parser.add_argument('-i', action="store", dest="iid", type=str, required=False, default="",
                        help="IID of the target individual. If unspecified, will use the prefix of the bam/mpileup file.")
    parser.add_argument('-t', action="store", dest="trim", type=int, required=False, default=0, 
                        help="trim certain number of bases from both ends.")
    parser.add_argument('-q', action="store", dest="q", type=int, required=False, default=30, 
                        help="Minimum mapping quality.")
    parser.add_argument('-Q', action="store", dest="Q", type=int, required=False, default=30, 
                        help="Minimum base quality.")
    parser.add_argument('--cleanup', action="store_true", dest="cleanup", required=False, 
                        help="whether to delete the intermediary hdf5 file. Default: False.")
    
    args = parser.parse_args()

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")
    from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom_BFGS

    conpop = []
    if args.conpop == None:
        conpop = ['CEU']
    else:
        conpop = [args.conpop]

    hapCon_chrom_BFGS(iid=args.iid, mpileup=args.mpileup, bam=args.bam, q=args.q, Q=args.Q,
    n_ref=2504, diploid_ref=False, exclude_pops=["AFR"], conPop=conpop, 
    h5_path1000g = args.ref, meta_path_ref = args.meta, folder_out="", 
    c=0.025, roh_jump=300, e_rate_ref=1e-3, logfile=False, output=False, cleanup=args.cleanup)
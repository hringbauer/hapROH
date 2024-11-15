import argparse
import h5py
from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom_BFGS
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Convert bam file to hdf5 format that stores readcount info at target sites.')
    parser.add_argument('-m', action="store", dest="mpileup", type=str, required=False,
                        help="path to samtools mpileup file")
    parser.add_argument('--bam', action="store", dest="bam", type=str, required=False,
                        help="path to bam file")
    parser.add_argument('--bamtable', action="store", dest="bamtable", type=str, required=False,
                        help="path to BamTable output")
    parser.add_argument('-r', action="store", dest="ref", type=str, required=True,
                        help="path to reference panel")
    parser.add_argument('--meta', action="store", dest="meta", type=str, required=True, 
                        default="/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv",
                        help="path to metadata of the reference panel")
    parser.add_argument('--con', action="store", dest="conpop", type=str, required=False, default="CEU",
                        help='source of contamination. Default is CEU. If you want to use allele frequencies calculated from multiple populations, use comma to delimit them. \
                            For example, --con CEU,TSI. There should be no whitespace surrouding the comma.')
    parser.add_argument('-c', action="store", dest="c", type=float, required=False, default=0.025,
                        help='Starting value for the optimization subroutine. Default is 0.025.')
    parser.add_argument('--exclude', action="store", dest='exHap', type=str, required=False, default="AFR",
                        help='Exclude certain haplotypes from the copying reference set. Default is AFR.\
                            If you do not want to exclude any haplotypes (this is especially important if you are working with aDNA from Africa),\
                                then you can simply provide a nonsense string that does not correspond to any populations in the 1000 Genome metadata. For example, --exclude 123.\
                                    If you want to exluce multiple subpopulations, use comma to separate strings. For example, --exclude CEU,TSI. There should be no whitespace surrounding the comma.')
    parser.add_argument('-i', action="store", dest="iid", type=str, required=False, default="",
                        help="IID of the target individual. If unspecified, will use the prefix of the bam/mpileup file.")
    parser.add_argument('--jump', action="store", dest="jump", type=float, required=False, default=300, 
                        help='Haplotype copying jump rate. Default is 300.')
    parser.add_argument('-e', action="store", dest="miscopy", type=float, required=False, default=1e-3, 
                        help='Haplotype copying error rate. This parameter models mutation/gene conversion/errors in the reference panel, etc. Default is 1e-3.')
    parser.add_argument('-q', action="store", dest="q", type=int, required=False, default=30, 
                        help="Minimum mapping quality. Only applicable when you use BAM file.")
    parser.add_argument('-Q', action="store", dest="Q", type=int, required=False, default=30, 
                        help="Minimum base quality. Only applicable when you use BAM file.")
    parser.add_argument('--lowmem', action="store_true", dest="lowmem", help="Use low memory mode.")
    parser.add_argument('-p', action="store", dest="prefix", type=str, required=False, default="hapCon")
    parser.add_argument('--log', action="store_true", dest="log",
                        help="Output a log file.")
    args = parser.parse_args()

    #sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")

    if args.conpop == 'OOA':
        conpop = ['EUR', 'EAS', 'SAS', 'AMR']
    elif ',' in args.conpop:
        conpop = [pop for pop in args.conpop.split(',')]
    else:
        conpop = [args.conpop]

    if ',' in args.exHap:
        exclude_pops = [pop for pop in args.exHap.split(',')]
    else:
        exclude_pops = [args.exHap]

    if args.lowmem:
        # check if the reference panel has binary-encoded genotypes
        with h5py.File(args.ref, 'r') as f:
            if not ('calldata' in f.keys() and 'GTbinary' in f['calldata'].keys()):
                raise RuntimeError(f'The reference panel {args.r} does not have binary-encoded genotypes. Please use a reference panel with the field calldata/binaryGT.')
        # check if the metadata of the ref panel has a column named "sex"
        meta = pd.read_csv(args.meta, sep="\t")
        if not "sex" in meta.columns:
            raise RuntimeError(f'The metadata of the reference panel {args.meta} does not have sex information')

    hapCon_chrom_BFGS(iid=args.iid, mpileup=args.mpileup, bam=args.bam, bamTable=args.bamtable, q=args.q, Q=args.Q,
            n_ref=2504, diploid_ref=True, exclude_pops=exclude_pops, conPop=conpop, 
            h5_path1000g = args.ref, meta_path_ref = args.meta, folder_out="", 
            c=args.c, roh_jump=args.jump, e_rate_ref=args.miscopy, lowmem=args.lowmem, 
            logfile=args.log, cleanup=False, prefix=args.prefix)
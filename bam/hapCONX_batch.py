import argparse
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert bam file to hdf5 format that stores readcount info at target sites.')
    parser.add_argument('--mlist', action="store", dest="mpileup", type=str, required=False,
                        help="path to a text file containing a list of mpileup files.")
    parser.add_argument('--blist', action="store", dest="bam", type=str, required=False,
                        help="path to a text file containing a list of BAM files.")
    parser.add_argument('-t', action="store", dest='t', type=int, required=False, default=1,
                        help="Number of processes to use. We recommend using multiple processes if your list contains multiple files.")
    parser.add_argument('-r', action="store", dest="ref", type=str, required=True,
                        help="path to reference panel")
    parser.add_argument('--meta', action="store", dest="meta", type=str, required=False, 
                        default="/mnt/archgen/users/yilei/Data/1000G/1000g1240khdf5/all1240/meta_df_all.csv",
                        help="path to metadata of the reference panel")
    parser.add_argument('--con', action="store", dest="conpop", type=str, required=False, default="CEU",
                        help='source of contamination. Default is CEU. If you want to use allele frequencies calculated from multiple populations, use comma to delimit them. \
                            For example, --con CEU, TSI.')
    parser.add_argument('-c', action="store", dest="c", type=float, required=False, default=0.025,
                        help='Starting value for the optimization subroutine. Default is 0.025.')
    parser.add_argument('--exclude', action="store", dest='exHap', type=str, required=False, default="AFR",
                        help='Exclude certain haplotypes from the copying reference set. Default is AFR.\
                            If you do not want to exclude any haplotypes (this is especially important if you are working with aDNA from Africa),\
                                then you can simply provide a nonsense string that does not correspond to any populations in the 1000 Genome metadata. For example, --exclude 123.\
                                    If you want to exluce multiple subpopulations, use comma to separate strings. For example, --exclude CEU, TSI.')
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
    parser.add_argument('-p', action="store", dest="prefix", type=str, required=False, default="hapCon")
    parser.add_argument('--log', action="store_true", dest="log",
                        help="Output a log file.")
    args = parser.parse_args()

    sys.path.insert(0, "/mnt/archgen/users/yilei/tools/hapROH/package")
    from hapsburg.PackagesSupport.hapsburg_run import hapCon_chrom_BFGS
    from hapsburg.PackagesSupport.parallel_runs.helper_functions import multi_run

    if args.mpileup and args.bam:
        print('please only provide one of BAM list or mpileup list...')
        sys.exit()

    sampleList = []
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

    if args.mpileup:
        with open(args.mpileup) as f:
            for line in f:
                iid, path2mpileup = line.strip().split()
                sampleList.append((iid, path2mpileup))
        prms = [[iid, path2mpileup, None, 30, 30, 2504, False, exclude_pops, conpop, 
                args.ref, args.meta, "", args.c, args.jump, args.miscopy,
                args.log, False, False, args.prefix] for iid, path2mpileup in sampleList]
    elif args.bam:
        with open(args.bam) as f:
            for line in f:
                iid, path2bam = line.strip().split()
                sampleList.append((iid, path2bam))
        prms = [[iid, None, path2bam, 30, 30, 2504, False, exclude_pops, conpop, 
                args.ref, args.meta, "", args.c, args.jump, args.miscopy,
                args.log, False, False, args.prefix] for iid, path2bam in sampleList]
    else:
        print('One of BAM or mpielup file list must be given.')
        sys.exit()
    
    print(f'running hapCon for {len(prms)} samples using {args.t} processes...')
    results = multi_run(hapCon_chrom_BFGS, prms, args.t)
    with open('hapCon.summary.tsv', 'w') as out:
        out.write(f'iid\tcontamX\tNumber_of_sites_covered\n')
        for sample, result in zip(sampleList, results):
            iid, _ = sample
            con_mle, lower, upper, numSitesCovered = result
            out.write(f'{iid}\t{round(con_mle,3)}({round(lower,3)} - {round(upper,3)})\t{numSitesCovered}\n')
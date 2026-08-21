"""Command-line interface for running hapROH's ROH calling functions.

Wraps `hapROH.run_new.callROH_chr`.

Examples:
    python hapROH.py -i sample.h5 -r ref_chr3.hdf5 -o results -c 3
    python hapROH.py -i sample.h5 -r ref_chr -o results -c 3
"""

import argparse

from hapROH.run_new import callROH_chr


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Call runs of homozygosity (ROH) using an HMM. ",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--path-sample",
        required=True,
        type=str,
        help="Path to the sample(s) file (Eigenstrat or HDF5) to call ROH on.",
    )
    parser.add_argument(
        "-r",
        "--path-ref",
        required=True,
        type=str,
        help="Path to the reference panel. If running on all chromosomes, suffix <chrom_number>.hdf5 will be appended",
    )
    parser.add_argument(
        "-o",
        "--folder-out",
        required=True,
        type=str,
        help="Base output directory. Results for each individual are written to folder_out/<iid>/chr<chrom>/.",
    )
    parser.add_argument(
        "-c",
        "--chrom",
        type=int,
        default=None,
        help="Chromosome number to process. If omitted, will run on all 22 autosomes.",
    )
    parser.add_argument(
        "--iid",
        dest="iids",
        nargs="+",
        default=None,
        metavar="IID",
        help="Individual ID(s) to call ROH for. If omitted, all individuals in the sample file are used.",
    )
    parser.add_argument(
        "--r-in",
        type=float,
        default=1,
        help="HMM transition rate into the ROH (homozygous) state. (default: %(default)s)",
    )
    parser.add_argument(
        "--r-out",
        type=float,
        default=20,
        help="HMM transition rate out of the ROH state. (default: %(default)s)",
    )
    parser.add_argument(
        "--r-jump",
        type=float,
        default=300,
        help="HMM jump rate between two distinct ROH states. (default: %(default)s)",
    )
    parser.add_argument(
        "--error-rate",
        type=float,
        default=0.01,
        help="Genotyping/sequencing error rate used by the emission model. (default: %(default)s)",
    )
    parser.add_argument(
        "--e-model",
        type=str,
        default="haploid",
        choices=["readcounts", "diploid_gt", "haploid"],
        help="Emission model to use. (default: %(default)s)",
    )
    parser.add_argument(
        "--downsampling",
        type=float,
        default=None,
        help="If provided, depth to downsample the sample data to before calling ROH. Only valid if input data is AD.",
    )
    parser.add_argument(
        "--cutoff-post",
        type=float,
        default=0.999,
        help="Cutoff used when calling ROH segments from the posterior probability. (default: %(default)s)",
    )
    parser.add_argument(
        "--snps-extend",
        type=int,
        default=0,
        help="Number of SNPs added to elongate ROH blocks (before merging). (default: %(default)s)",
    )
    parser.add_argument(
        "--max-gap",
        type=float,
        default=0.005,
        help="Maximum gap (in Morgans) between two adjacent ROH for them to be merged. (default: %(default)s)",
    )
    parser.add_argument(
        "--min-len1",
        type=float,
        default=0.04,
        help="Minimum length (in Morgans) required for two adjacent ROH to be merged. (default: %(default)s)",
    )
    parser.add_argument(
        "--min-len2",
        type=float,
        default=0.02,
        help="Minimum length (in Morgans) required for two adjacent ROH to be merged. (default: %(default)s)",
    )
    parser.add_argument(
        "--min-len-final",
        type=float,
        default=0.04,
        help="Minimum length (in Morgans) for segments to appear in the final dataset. (default: %(default)s)",
    )
    parser.add_argument(
        "--logfile",
        type=str,
        default=None,
        help="Path to a file to write log output to. If omitted, logs go to the default stream handler.",
    )
    parser.add_argument(
        "--loglevel",
        type=int,
        default=1,
        help="Verbosity level for the hapROH logger (0=WARNING, 1=INFO, 2 or higher=DEBUG). (default: %(default)s)",
    )
    parser.add_argument(
        "--backend",
        type=str,
        default="cython",
        choices=["python", "numba", "cython"],
        help="Which backend to the for the main computation algorithm.",
    )

    args = parser.parse_args()
    kwargs = vars(args)
    chrom_arg = kwargs.pop("chrom")

    if chrom_arg is not None:
        callROH_chr(chrom=chrom_arg, **kwargs)
    else:
        base_path_ref = kwargs["path_ref"]
        for chrom in range(1, 23):
            kwargs["path_ref"] = f"{base_path_ref}{chrom}.hdf5"
            callROH_chr(chrom=chrom, **kwargs)


if __name__ == "__main__":
    main()

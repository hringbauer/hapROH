#!/usr/bin/env python
"""Command-line interface to prepare a BAM file into a hdf5 file using samtools mpileup

Wraps `hapROH.utils.IO.bam2hdf5`.

Examples:
    python bam_to_hdf5.py --path-bam sample.bam --path-ref ref_chr3.hdf5 \\
        --chrom 3 --path-out sample_chr3.hdf5

    python bam_to_hdf5.py --path-bam sample.bam --path-ref ref_chr \\
        --path-out sample_chr
"""
import argparse

from hapROH.utils.IO import bam2hdf5

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert a bam into a hdf5 file using samtools mpileup",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--path-bam", required=True, type=str,
                         help="Path to the bam file.")
    parser.add_argument("-r", "--path-ref", required=True, type=str,
                         help="Path to the reference panel, defining which positions are extracted. If running on all chromosomes, suffix <chrom_number>.hdf5 will be appended")
    parser.add_argument("-o", "--path-out", required=True, type=str,
                         help="Path out. If running on all chromosomes, suffix <chrom_number>.hdf5 will be appended")
    parser.add_argument("--sample-name", required=True, type=str,
                         help="Name to use in the hdf5.")
    parser.add_argument("-c", "--chrom", type=int, default=None,
                         help="Chromosome number to process. If omitted, will run on all 22 autosomes.")
    parser.add_argument("-Q", "--min-BQ", dest="min_base_qual", type=int, default=30,
                         help="Minimum base quality for a base to be considered.")
    parser.add_argument("-q", "--min-MQ", dest="min_map_qual", type=int, default=None,
                         help="Minimum mapping quality for an alignment to be used")
    parser.add_argument("--overwrite", action="store_true",
                         help="Overwrite existing bam files.")
    parser.add_argument("--samtools", type=str, default="samtools",
                         help="Name of the command or executable path for samtools")

    args = parser.parse_args()
    kwargs = vars(args)
    chrom_arg = kwargs.pop("chrom")

    # Verify the samtools command/executable actually works
    samtools_cmd = kwargs["samtools"]
    try:
        result = subprocess.run(
            [samtools_cmd, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except (OSError, FileNotFoundError) as e:
        raise RuntimeError(
            f"Could not run samtools using command/path '{samtools_cmd}': {e}. "
            "Please check that samtools is installed and that --samtools points "
            "to a valid command name or executable path."
        ) from e

    if chrom_arg is not None:
        bam2hdf5(chrom=chrom_arg, **kwargs)
    else:
        base_path_ref = kwargs["path_ref"]
        base_path_out = kwargs["path_out"]
        for chrom in range(1, 23):
            kwargs["path_ref"] = f"{base_path_ref}{chrom}.hdf5"
            kwargs["path_out"] = f"{base_path_out}{chrom}.hdf5"
            bam2hdf5(chrom=chrom, **kwargs)

if __name__ == "__main__":
    main()

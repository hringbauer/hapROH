#!/usr/bin/env python
"""Command-line interface to prepare a BAM file into a hdf5 file using samtools mpileup

Wraps `hapROH.utils.IO.bam2hdf5`.

Examples:
    python bam_to_hdf5.py -i sample.bam -o sample_chr3.hdf5 -r ref_chr3.hdf5 -c 3 -s my_sample
    python bam_to_hdf5.py -i sample.bam -o sample_chr -r ref_chr -s my_sample
"""
import argparse, subprocess

from hapROH.utils.IO import bam2hdf5

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert a bam into a hdf5 file using samtools mpileup",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--path-bam", required=True, type=str,
                         help="Path to the input bam file.")
    parser.add_argument("-r", "--path-ref", dest="path_refHDF5", required=True, type=str,
                         help="Path to the reference HDF5, defining which positions are extracted. If running on all chromosomes, suffix <chrom_number>.hdf5 will be appended")
    parser.add_argument("-o", "--path-out", dest="path_outHDF5", required=True, type=str,
                         help="Path to the output HDF5. If running on all chromosomes, suffix <chrom_number>.hdf5 will be appended")
    parser.add_argument("-s", "--sample-name", required=True, type=str,
                         help="Name to use in the hdf5.")
    parser.add_argument("-c", "--chrom", type=int, default=None,
                         help="Chromosome number to process. If omitted, will run on all 22 autosomes.")
    parser.add_argument("-Q", "--min-BQ", dest="min_base_qual", type=int, default=30,
                         help="Minimum base quality for a base to be considered.")
    parser.add_argument("-q", "--min-MQ", dest="min_map_qual", type=int, default=None,
                         help="Minimum mapping quality for an alignment to be used")
    parser.add_argument("--overwrite", action="store_true",
                         help="Overwrite existing bam files.")
    parser.add_argument("--samtools", dest="path_samtools", type=str, default="samtools",
                         help="Name of the command or executable path for samtools")

    args = parser.parse_args()
    kwargs = vars(args)
    chrom_arg = kwargs.pop("chrom")

    # Verify that the samtools command/executable actually works
    samtools_cmd = kwargs["path_samtools"]
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
        base_path_ref = kwargs["path_refHDF5"]
        base_path_out = kwargs["path_outHDF5"]
        for chrom in range(1, 23):
            kwargs["path_refHDF5"] = f"{base_path_ref}{chrom}.hdf5"
            kwargs["path_outHDF5"] = f"{base_path_out}{chrom}.hdf5"
            bam2hdf5(chrom=chrom, **kwargs)

if __name__ == "__main__":
    main()

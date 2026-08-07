### Python Functions to prepare and manipulate input files
### 2026

### Imports
import os, subprocess, tempfile, logging
import numpy as np
import pandas as pd
import h5py
import re

logger = logging.getLogger(__name__)

#########################
# Utility functions to convert bam to hdf5
# by Florence P., developed on Leipzig MPI Server, July 2026

def get_snp_from_h5(path_h5: str) -> pd.DataFrame:
    """Extract biallelic SNP positions and alleles from a reference HDF5 panel.

    Args:
        path_h5: Path to the reference HDF5 file.

    Returns:
        A DataFrame containing the biallelic SNP sites. Its columns are ``pos``,
        ``map``, ``ref``, ``alt``, and, when available in the input file,
        ``chrom``.
    """
    with h5py.File(path_h5, "r") as h5_ref:
        variants = h5_ref["variants"]
        pos = np.array(variants["POS"]).astype(int)
        map = np.array(variants["MAP"]).astype(float)
        ref = np.array(variants["REF"]).astype("U1")
        alt = np.array(variants["ALT"]).astype("U1")
        chrom = (
            np.array(variants["CHROM"]).astype(int)
            if "CHROM" in variants
            else None
        )

    # filter only biallelic variants and reshape alt into a 1d-array if needed
    mask_biallelic = np.ones(len(alt), dtype=bool)
    if alt.ndim == 2:
        mask_biallelic = np.all(alt[:, 1:] == "", axis=1)
        alt = alt[:, 0]

    # filter only SNP
    bases = np.array(["A", "T", "G", "C"])
    mask_snp = np.isin(ref, bases) & np.isin(alt, bases)

    mask = mask_biallelic & mask_snp
    print(f"\tKept {np.sum(mask)}/{len(mask)} biallelic SNP sites")

    df_snp = pd.DataFrame({
        "pos": pos,
        "map": map,
        "ref": ref,
        "alt": alt,
    })

    if chrom is not None:
        df_snp["chrom"] = chrom

    return df_snp[mask]

def bam2pileup(path_bam:str, path_refHDF5:str, chrom:int|None=None, min_base_qual:int=30, min_map_qual:int=30, path_samtools="samtools") -> pd.DataFrame:
    """Count reference and alternate alleles in a BAM file at reference-panel SNPs.

    Args:
        path_bam: Path to the input BAM file.
        path_refHDF5: Path to the reference HDF5 panel that defines SNP sites.
        chrom: Chromosome to process. Required when the reference
            HDF5 file does not contain chromosome information.
        min_base_qual: Minimum base quality passed to ``samtools mpileup``.
        min_map_qual: Minimum mapping quality passed to ``samtools mpileup``.
        path_samtools: Path or command name for the ``samtools`` executable.

    Returns:
        A DataFrame with the columns ``chrom``, ``pos``, ``map``, ``ref``,
        ``alt``, ``ref_count``, and ``alt_count``.
    """
    bases = np.array(['A', 'T', 'G', 'C'])

    ### Load SNPs from reference hdf5
    print(f"Extracting SNP positions from hdf5")
    df_h5 = get_snp_from_h5(path_refHDF5)
    if not "chrom" in df_h5.columns:
        if chrom is None:
            raise ValueError("The reference hdf5 does not contain a field 'chrom'. Please specify it as an argument")
        else:
            df_h5["chrom"] = chrom
    elif chrom is not None:
        mask = df_h5["chrom"] == chrom
        if mask.sum() == 0:
            raise ValueError(f"Contig {chrom} not found in hdf5 (chromosomes in hdf5: {df_h5["chrom"].unique()})")
        df_h5 = df_h5[mask]

    ### Run pileup on BAM file with temporary bed file
    print(f"Running mpileup on BAM file")
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create temp bedfile to pileup ony desired positions:
        df_bed = pd.DataFrame({"chrom":df_h5["chrom"], "start":df_h5["pos"]-1, "end":df_h5["pos"]})
        path_bed = os.path.join(tmp_dir, "positions.bed")
        df_bed.to_csv(path_bed, header=False, index=False, sep="\t")

        # launch samtools pileup as a subprocess
        command_pileup = [
            path_samtools, "mpileup", "--no-BAQ", "-Q", str(min_base_qual), "-q", str(min_map_qual),
            path_bam, "--positions", path_bed,
        ]
        proc = subprocess.Popen(
            command_pileup,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        # pipe pileup output into pandas 
        df_pileup = pd.read_csv(
            proc.stdout,
            sep="\t",
            header=None,
            names=["chrom", "pos", "ref", "depth", "read_bases", "base_qualities"],
        )
        returncode = proc.wait()
        stderr_output = proc.stderr.read()

        if returncode != 0:
            raise RuntimeError("Pileup command failed with return code: {returncode}", stderr_output)

    ### Process output from pileup
    mask = df_pileup["depth"] > 0
    df_pileup = df_pileup[mask]
    print(f"\t {sum(mask)}/{len(mask)} positions covered, at mean depth {np.mean(df_pileup["depth"])}")

    # add SNP data (map, ref, alt) from hdf5
    df_pileup = df_pileup.drop("ref", axis=1).merge(df_h5, on=["chrom", "pos"])

    # ignore positions with deletions (*#) and skips (<>) or adjacent to indels (+-), 
    characters = "+-<>*#"
    re_exp = "|".join([re.escape(c)for c in characters])
    mask = df_pileup["read_bases"].str.contains(re_exp, regex=True)
    df_pileup = df_pileup[~mask]
    print(f"\t Ignored {sum(mask)}/{len(mask)} positions containing indels")

    # remove ^X (read start + mapping quality char) and $ (read end)
    df_pileup["clean"] = (
        df_pileup["read_bases"]
        .str.replace(r"\^.", "", regex=True)  # start-of-read marker '^' plus the following mapqual char
        .str.replace(r"\$", "", regex=True)   # end-of-read marker '$'
        .str.upper()
    )

    # get matrix of read count
    count_matrix = []
    for base in bases:
        count_matrix.append(df_pileup["clean"].str.count(base))
    count_matrix = np.array(count_matrix)

    # get idx of column containing ref and alt count for each row
    base2idx = dict([(c,i) for i, c in enumerate(bases)])
    ref_idx = df_pileup["ref"].str.upper().map(base2idx).to_numpy()
    alt_idx = df_pileup["alt"].str.upper().map(base2idx).to_numpy()

    rows = np.arange(len(df_pileup))
    df_pileup["ref_count"] = count_matrix[ref_idx, rows]
    df_pileup["alt_count"] = count_matrix[alt_idx, rows]

    ### Check if bases correspond to expected ref/alt ones
    total_bases = sum(df_pileup["depth"])
    correct_bases = sum(df_pileup["ref_count"] + df_pileup["alt_count"])
    print(f"\t Found {total_bases-correct_bases}/{total_bases}={100*(total_bases-correct_bases)/total_bases:.3}% bases different from expected ref/alt")  # =2/3 sequencing error rate ?
    # remove positions where the base does not correspond to the expected ones
    mask = df_pileup["depth"] > df_pileup["ref_count"] + df_pileup["alt_count"]
    df_pileup = df_pileup[~mask]
    print(f"\t Ignored {sum(mask)}/{len(mask)} positions with bases different from expected ref/alt")

    # cleanup intermediate columns
    df_pileup = df_pileup.drop(columns=["read_bases", "base_qualities", "clean", "depth"])

    return df_pileup

def pileup2hdf5(df_pileup:pd.DataFrame, path_outHDF5:str, sample_name:str, overwrite:bool=False) -> None:
    """Write pileup allele counts to an HDF5 file.

    Args:
        df_pileup: DataFrame produced by :func:`bam2pileup`.
        path_outHDF5: Path for the output HDF5 file.
        sample_name: Sample identifier to store in the output file.
        overwrite: Whether to overwrite or fail if the file already exists.

    Returns:
        None.
    """
    os.makedirs(os.path.dirname(path_outHDF5), exist_ok=True)
    with h5py.File(path_outHDF5, 'w' if overwrite else 'w-') as f_out:
        f_out.create_dataset("variants/CHROM", data=df_pileup["chrom"])
        f_out.create_dataset("variants/POS", data=df_pileup["pos"])
        f_out.create_dataset("variants/MAP", data=df_pileup["map"])
        f_out.create_dataset("variants/REF", data=df_pileup["ref"].astype('S1'))
        f_out.create_dataset("variants/ALT", data=df_pileup["alt"].astype('S1'))
        f_out.create_dataset("calldata/AD", data=df_pileup[["ref_count", "alt_count"]].to_numpy()[:, np.newaxis, :]) # shape (nb_snp, nb_samples=1, 2)
        f_out.create_dataset("samples", data=np.array([sample_name]).astype('S50'))

def bam2hdf5(path_bam:str, path_refHDF5:str, path_outHDF5:str, sample_name:str, chrom:int|None=None, overwrite:bool=False, min_base_qual:int=30, min_map_qual:int=30, path_samtools="samtools") -> None:
    """Convert a BAM file to an HDF5 file using a reference SNP panel. Final HDF5 file will have following field: samples; variants/POS,MAP,REF,ALT,[CHROM]; calldata/AD

    Args:
        path_bam: Path to the input BAM file.
        path_refHDF5: Path to the reference HDF5 panel that defines SNP sites.
        path_outHDF5: Path for the output HDF5 file.
        sample_name: Sample identifier to store in the output file.
        chrom: Chromosome or contig to process. Required when the reference
            HDF5 file does not contain chromosome information.
        overwrite: Whether to overwrite or fail if the file already exists.
        min_base_qual: Minimum base quality passed to ``samtools mpileup``.
        min_map_qual: Minimum mapping quality passed to ``samtools mpileup``.
        path_samtools: Path or command name for the ``samtools`` executable.

    Returns:
        None.
    """
    if not overwrite and os.path.isfile(path_outHDF5):
        print(f"File {path_outHDF5} exists already and overwrite=False. Nothing happend.")
        return
    df_pileup = bam2pileup(path_bam, path_refHDF5, chrom, min_base_qual, min_map_qual, path_samtools)
    pileup2hdf5(df_pileup, path_outHDF5, sample_name, overwrite)

def bam2hdf5s(path_bam:str, prefix_refHDF5:str, dir_outHDF5:str, sample_name:str, overwrite:bool=False, min_base_qual:int=30, min_map_qual:int=30, path_samtools="samtools") -> None:
    """Wrapper to call ``bam2hdf5`` on each chromosome from 1 to 22.

    Args:
        path_bam: Path to the input BAM file.
        prefix_refHDF5: Prefix of reference HDF5 paths; chromosome suffixes
            from ``1.hdf5`` through ``22.hdf5`` are appended.
        prefix_outHDF5: Output directory, in which HDF5 files of format $iid.chr$ch.hdf5 will be created
        sample_name: Sample identifier to store in each output file.
        overwrite: Whether to overwrite or fail if the files already exist.
        min_base_qual: Minimum base quality passed to ``samtools mpileup``.
        min_map_qual: Minimum mapping quality passed to ``samtools mpileup``.
        path_samtools: Path or command name for the ``samtools`` executable.

    Returns:
        None.
    """
    prefix_outHDF5 = os.path.join(dir_outHDF5, sample_name+".chr")
    for chrom in range(1, 23):
        suffix = str(chrom) + ".hdf5"
        bam2hdf5(path_bam, prefix_refHDF5+suffix, prefix_outHDF5+suffix, sample_name, chrom, overwrite, min_base_qual, min_map_qual, path_samtools)

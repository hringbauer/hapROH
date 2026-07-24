### Python Functions to prepare and manipulate input files
### 2026

### Imports
import os, warnings, subprocess, tempfile
import numpy as np
import pandas as pd
import h5py



### Functions Utils1 and Utils2 by Florence P., developed on Leipzig MPI Server

### Utils 1: Create the .snp from the reference pannel

def get_snp_from_h5(path_h5:str, path_snp:str="", path_bed:str="", overwrite:bool=False):
    """Extract the position and information of the SNPs present in the reference hdf5 panel, 
    and saves them into eigenstrat (.snp) or bed (.bed) format"""
    for path in [path_snp, path_bed]:
        if os.path.exists(path):
            if not overwrite:
                warnings.warn(f"Warning: output file {path} exists already. Set overwrite to True to continue")
                return
            else:
                os.remove(path)
    for chrom in range(1, 23):
        path_h5_ch = path_h5 + str(chrom) + ".hdf5"

        with h5py.File(path_h5_ch, "r") as h5_ref:
            posBP = np.array(h5_ref["variants/POS"]).astype(str)
            posM = np.array(h5_ref["variants/MAP"]).astype(str)
            ref = np.array(h5_ref["variants/REF"]).astype(str)
            alt = np.array(h5_ref["variants/ALT"]).astype(str)

        df_snp = pd.DataFrame({
            "posM": posM,
            "posBP": posBP,
            "ref": ref,
            "alt": alt})
        
        df_snp["chr"] = str(chrom)
        df_snp["snp_id"] = df_snp["chr"] + ":" + df_snp["posBP"]
        df_snp = df_snp[["snp_id", "chr", "posM", "posBP", "ref", "alt"]]

        if path_snp != "":
            df_snp.to_csv(path_snp, mode="a", index=False, header=False, sep='\t')

        if path_bed != "":
            df_bed = df_snp[["chr"]]
            df_bed["start"] = df_snp["posBP"].astype(int) - 1
            df_bed["end"] = df_snp["posBP"]
            df_bed.to_csv(path_bed, mode="a", index=False, header=False, sep='\t')

### Utils 2: Create the eigenstrat from the bam file

def bam2eigenstrat(path_bam: str, sample_name: str, path_out:str="", 
                   path_h5: str = "", path_snp: str = "", q=30, Q=30,
                   pileupcaller="/home/harald_ringbauer/.local/bin/pileupCaller"):
    """Creates a pseudo-haploid eigenstrat file from a bam file. SNP positions are provided
    either by the reference h5 panel or by a .snp file (eigenstrat format).

    Args:
        path_bam: path to the BAM file to convert
        sample_name: name of the sample, will be used in .ind file and to name the output files
        path_out: path to the outputs files, without file extension (.ind, .snp, .geno)
        path_h5: path to the reference panel in hdf5 format, up to XX.h5 (with XX the chromosome
            number). Ignored if path_snp is provided
        path_snp: path to the snp to use, in eigenstrat (.snp) format
        q,Q: Parameters for samtools pileup (q minimum mapping quality,Q minimum base quality
        pileupCaller: Path to the pileupcaller binary. If empty, use the command pileupCaller
    """
    ### Use default pileupcaller command if none given
    if len(pileupcaller)==0:
        pileupcaller="pileupCaller"
        
    try:
        subprocess.run(pileupcaller, capture_output=True)
    except FileNotFoundError:
        warnings.warn("pileupCaller is missing")

    if path_snp == "" and path_h5 == "":
        raise ValueError("You must provide either path_snp or path_h5")

    if path_out == "":
        path_out = sample_name

    ### Run pileup and then pileupcaller
    with tempfile.TemporaryDirectory() as tmp_dir:
    #tmp_dir = "/mnt/archgen/users/hringbauer/data/med_jews/"
    
        if path_snp == "":
            path_snp = os.path.join(tmp_dir, "positions.snp")   # around 25s
            get_snp_from_h5(path_h5, path_snp)
    
        path_pos = os.path.join(tmp_dir, "positions.pos")
        df_snp = pd.read_csv(path_snp, header=None, sep="\t")
        df_pos = df_snp[[1, 3]]
        df_pos.to_csv(path_pos, header=False, index=False, sep="\t")
    
        command_pileup = [                                      # around 8s
            "samtools", "mpileup", "--no-BAQ", f"-q30", f"-Q30",
            path_bam, "--positions", path_pos]
            
        command_caller = [                                      # around 150s
            pileupcaller, "--randomHaploid",
            "--sampleNames", sample_name,
            "--snpFile", path_snp,
            "-e", path_out]
    
        with subprocess.Popen(
            command_pileup, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ) as pileup_proc:
            with subprocess.Popen(
                command_caller,
                stdin=pileup_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            ) as caller_proc:
                pileup_proc.stdout.close()
                caller_out, caller_err = caller_proc.communicate()
    
            _, pileup_err = pileup_proc.communicate()
    
        if pileup_proc.returncode != 0:
            raise RuntimeError(f"samtools mpileup failed:\n{pileup_err}\n{caller_err}")
        if caller_proc.returncode != 0:
            raise RuntimeError(f"pileupCaller failed:\n{caller_err}")

    print(f"Done !\nOutput written to {sample_name}.snp/.ind/.geno")
    return
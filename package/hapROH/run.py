### Key Python Functions to run hapROH
### Created by Harald Ringbauer hapROH v0.67
### To bundle hapROH run imports more efficiently and clear for the user
### Wraps the "original" functions in hapsburg.PackagesSupport
### and makes them accessible via shorter import statements.

from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind_lowmem
from hapsburg.PackagesSupport.pp_individual_roh_csvs import pp_individual_roh

def callROH_ind(iid, chs=range(1,23),
              path_targets_prefix="", path_targets="",
              h5_path1000g=None, meta_path_ref=None, folder_out=None, prefix_out="",
              e_model="haploid", p_model="Eigenstrat", post_model="Standard", low_mem=True,
              processes=1, delete=False, verbose=True, save=True, save_fp=False, 
              n_ref=2504, diploid_ref=True, exclude_pops=[], readcounts=True, random_allele=True, downsample=False,
              c=0.0, conPop=["CEU"], roh_in=1, roh_out=20, roh_jump=300, e_rate=0.01, e_rate_ref=0.00, 
              cutoff_post = 0.999, max_gap=0.005, roh_min_l_initial = 0.02, roh_min_l_final = 0.04,
              min_len1 = 0.02, min_len2 = 0.04, logfile=True, combine=True, file_result="_roh_full.csv"):
    """Main hapROH function to call autosomal ROH. 
    Analyze a full single individual in a parallelized fashion. Run multiple chromosome analyses in parallel.
    Bring together the result ROH tables from each chromosome into a genome-wide summary ROH table.
    This function wraps hapsb_chrom. The default Parameters are finetuned for pseudo-haploid 1240k aDNA data.

    Parameters
    ----------
    iid: str
        IID of the Target Individual, as found in Eigenstrat.
    chs: list
        Set of chromosomes to call ROH on
    path_targets_prefix:
        A directory containing an HDF5 file for each chromosome. The file name should follow $iid.chr$ch.hdf5.
    path_targets: str
        Path of the target files. You need only specify one of path_targets_prefix or path_targets.
    h5_path1000g: str
        Path of the reference genotypes
    meta_path_ref: str 
        Path of the meta file for the references
    folder_out: str
        Path of the basis folder for output
    prefix_out: str
        Path to insert in output string, e.g. test/ [str]
    e_model: str
        Specify Emission model to use, should be one of haploid/diploid_gt/readcount
    p_model: str
        Specify Preprocessing model, should be one of EigenstratPacked/EigenstratUnpacked/MosaicHDF5
    post_model: str
        Specify post-process model, should be one of Standard/MMR (experimental)
    low_mem:
        Whether to use the low memory version of hapROH (default since v0.67).
    processes: int
        How many Processes to use
    delete: bool
        Whether to delete raw posterior per locus
    verbose: bool
        Whether to print extensive output
    save: bool
        Whether to save the inferred ROH
    save_fp: bool
        Whether to save the full posterior matrix
    n_ref: int
        Number of (diploid) reference Individuals to use
    diploid: bool
        Whether the reference panel is diploid or not (e.g., for an autosome, True, and for the male X chromosome, False). 
    exclude_pops: list of str
        Which populations to exclude from the reference panel
    readcounts: bool
        Whether to load readcount data
    random_allele: bool
        Whether to pick a random one of the two target alleles per locus
    downsample:
        If not false (i.e., float), downsample readcounts to this target average coverage
    c: float
        Contamination rate. This is only applicable if the emission model is readcount_contam.
    conPop: list of str
        Ancestry of contamination source. Only applicable if the emission model is readcount_contam.
    roh_in: float
        Parater to jump into ROH state (per Morgan)
    roh_out: float
        Parameter to jump out of ROH state (per Morgan)
    roh_jump: float
        Parameter to jump (per Morgan)
    e_rate: float
        Sequencing error rate.
    e_rate_ref: float
        Haplotype miscopying rate.
    cutoff_post: float
        Posterior cutoff for ROH calling
    max_gap: float
        Maximum gap to merge two adjacent short ROH blocks (in Morgan)
    roh_min_l_initial: float
        Minimum length of ROH blocks to use before merging adjacent ones (in Morgan)
    roh_min_l_final: float
        Minimum length of ROH blcoks to output after merging (in Morgan)
    min_len1: float
        Minimum length of the shorter candidate block in two adjacent blocks that can be merged (in Morgan)
    min_len2: float
        Minimum length of the longer candidate block in two adjacent blocks that can be merged (in Morgan)
    logfile: bool
        Whether to use the logfile (or print to sys out)
    combine: bool 
        Whether to combine the output of all chromosomes
    file_result: str
        Appendix to individual results

    Return: If combine is true, return a pandas dataframe that contains information of all detected ROH blocks. Otherwise, nothing is returned.
    """ 

    ### Switch to the new (since v0.67, more memory-efficient and faster) or the old hapROH algorithm
    if low_mem:
        hapsb_ind_lowmem(iid=iid, chs=chs, path_targets_prefix=path_targets_prefix, path_targets=path_targets,
                  h5_path1000g=h5_path1000g, meta_path_ref=meta_path_ref, folder_out=folder_out, prefix_out=prefix_out,
                  e_model=e_model, p_model=p_model, post_model=post_model,
                  processes=processes, delete=delete, output=verbose, save=save, save_fp=save_fp, 
                  n_ref=n_ref, diploid_ref=diploid_ref, exclude_pops=exclude_pops, readcounts=readcounts, random_allele=random_allele, downsample=downsample,
                  c=c, conPop=conPop, roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump, e_rate=e_rate, e_rate_ref=e_rate_ref, 
                  cutoff_post = cutoff_post, max_gap=max_gap, roh_min_l_initial = roh_min_l_initial, roh_min_l_final = roh_min_l_final,
                  min_len1 = min_len1, min_len2 = min_len2, verbose=verbose, logfile=logfile, combine=combine, file_result=file_result)

    else:
        hapsb_ind(iid=iid, chs=chs,
                  path_targets_prefix=path_targets_prefix, path_targets=path_targets,
                  h5_path1000g=h5_path1000g, meta_path_ref=meta_path_ref, folder_out=folder_out, prefix_out=prefix_out,
                  e_model=e_model, p_model=p_model, post_model=post_model,
                  processes=processes, delete=delete, output=verbose, save=save, save_fp=save_fp, 
                  n_ref=n_ref, diploid_ref=diploid_ref, exclude_pops=exclude_pops, readcounts=readcounts, random_allele=random_allele, downsample=downsample,
                  c=c, conPop=conPop, roh_in=roh_in, roh_out=roh_out, roh_jump=roh_jump, e_rate=e_rate, e_rate_ref=e_rate_ref,
                  cutoff_post = cutoff_post, max_gap=max_gap, roh_min_l_initial = roh_min_l_initial, roh_min_l_final = roh_min_l_final,
                  min_len1 = min_len1, min_len2 = min_len2, logfile=logfile, combine=combine, file_result=file_result)    


def ppROH_inds(iids, base_folder="./Output/", suffix='_roh_full.csv', 
               save_path="", meta_path="", 
               min_cm=[4,8,12,20], snp_cm=50, gap=0.5, min_len1=2.0, min_len2=4.0, 
               verbose=True, meta_info=True):
    """Post-process Individual ROH .csv files. Combines them into a single summary ROH.csv file, saved in save_path.
    Use Individuals' IDs, create paths, and run the combining.
    iids: List of target Individuals
    base_folder: Folder where to find individual results .csvs
    min_cm: Minimum post-processed Length of ROH blocks. Array (to have multiple possible values)
    snp_cm: Minimum Number of SNPs per cM
    gap: Maximum length of gaps to merge
    output: Whether to plot output per Individual.
    meta_info: Whether to merge in Meta-Info from the original Meta File
    save_path: If given, save the resulting dataframe there
    meta_path: [Optional] Path to Meta data table. If provided, merge that table in.
    Should contain clst column - this sets the population labels (also later used for plotting)
    min_len1: Minimum Length of shorter block to merge [cM]
    min_len2: Maximum Length of longer block to merge [cM]"""

    df = pp_individual_roh(iids, meta_path=meta_path, base_folder=base_folder, 
                      suffix=suffix, save_path=save_path, min_cm=min_cm, snp_cm=snp_cm, 
                      gap=gap, min_len1=min_len1, min_len2=min_len2, output=verbose, meta_info=meta_info)

    return df
    






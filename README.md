# hapROH
Software to call ROH from ancient and present-day DNA using reference haplotypes.
Author: Harald Ringbauer, September 2020
Code released under GNU General Public License v3.0

Development Code behind python package hapROH.
See https://pypi.org/project/hapROH/ for detail about the current version and installation.

This is the git repository for development, as well as code used in the hapROH publication
(Preprint at https://doi.org/10.1101/2020.05.31.126912)

This code is experimental, and desgined to run on the University of Chicago midway cluster.

The code behind the hapROH package can be found in `./package/`

I do not share the data here, so most of the code in these notebooks will not run. The code is shared for transparency reason. If you need some of the data, feel free to reach out directly to me.

# Internal notes about folder structure:

## Data Preparation:
For 1240k reference data:
There is a notebook `prepare_1000genomes.ipynb`:
It downloads 1000 Genomes to local machine and prepares
downsampled and clean 1240k data (subsampled to biallelic SNPs), and saves it to hdf5

## Running the HMM

### Running Modes:
cython=
0: Python Fwd/Bwd Algorithm
1: Full Cython Algorithm
2: Optimized Cython Algorithm (linear Nr. references, full Linkage Map Model)

# Folders:
## Data
./Data Folder is git-ignored. That's where the data goes into.

## /1000 Genome Folder Structure

### /Individuals
Contains one column lists of individuals to use

### /Markers
Contains 1240k Eigenstrat .snp, that specifies genomic position and map position of every SNP

### /Autosome VCF
Contains the input VCFs, downloaded from 1000 Genomes

### Autosome VCF/Subset
Contains also the processed VCFs, downsampled to 

## HDF5
### \FULLHDF5
Contains the full HDF5s (used temporarily to create downsampled HDF5s)

### \1240kHDF5
Contains the final product. Downsampled hdf5s to 1240k biallelic SNPs

## /ReichLabEigenstrat
Contains the Eigenstrat from the Reich website. Use script in ./Notebooks/PrepareData/prepare_Reich_Eigenstrat to download and unzip it



# How to run Notebooks on cluster:
sbatch ./jupyter.sbatch

and then collect the URL from the .err file of the job. Afterwards: Do not forget to scancel the job otherwise it runs until the time limit is hit.

This allows one to run memory or cpu-intense tasks via a notebook on the cluster.

For debugging/info: Show job account information for a specific job:
sacct -j jobid --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss,MaxVMSize,nnodes,ncpus,nodelist



## Notebooks to run for creation of data used in paper
Several Notebooks to run HAPSBURG in parallel on simulated as well as empirical data are found in
`Notebooks/ParallelRuns/`

these run Hapsburg in a parallelized fashion, assuming that the notebook is run on a cluster with enough nodes (which can be set). The general strategy is to produce lists of parameter files (with one list of parameters per run, e.g. which individual and chromosome), which are then spread out by 'pool.starmap' to nodes via a multirun function.

## Also: HAPSBURG  is run on cluster with slurm sequentially
Code can be found in ./scripts/cluster_runs/

Usually they are written as array jobs:
Each folder contains a python file that runs a single individual, taking an integer as input. These integers are mapped to individuals using the accoring meta files (somewhere in ./data/)

The second file in each folder is `slurm_batch.sh` which is the shell submission script. Modify the .sh scripts to specify what batch to run (i.e. how many and which individuals to submit)

They are submitted with:
`sbatch slurm_batch.sh`

Usually Analysis per chromosomes are run, and results saved into folders: basefolder/iid/chr0/

There is a helper function `combine_individual_data` from `./packagesSupport/parallel_runs/helper_functions` that combines these results per Chromosome into a .csv per Individual, with a option to delete the Chromosome Result folder (to save space). For detailled downstream analysis (e.g. of the Posterior), the `delete` argument has to be False to keep the posterior file. The standard output ist `base_folder/IID_roch_full.csv`

## Further Processing per Individual:
These jobs produce summary files for each individual (ususally with `iid_roh_full.csv`) with fields
(Start,End,StartM,EndM,length,lengthM,iid,ch)

These are then further post-process with functions in `./packagesSupport/pp_individual_roh_csvs.py`.

The ususual structure of doing this is:

1) Produce list of paths to .csvs
2) Iterate over these csvs, extract and postprocess (e.g. gap-merging) these
3) These are then merged into one master .csv file with the fields, saved as `combinedroh.csv` with columns
(iid, pop, max_roh, sum_roh, n_roh, lat, lon, age), where some summary details can be merged back in from the master meta dataframe. The function to do so can be found in pp_individual_roh_csvs.py in PackagesSupport. This function is typically called in the ParallelRun notebooks, towards the end

One example of this workflow: `notebooks/PrepareData/roh_process.ipynb`

### Combining Different Dataset
These combined csvs are then combined again in 'combine_roh_results.ipynb' into one giant master dataframe. In this function, also "region" columns are defined, which are the basis for downstream processing.

### Example: Antonio 2019 Individuals
There is a notebook in `Notebooks/PrepareData` to prepare the Meta, as well as the hdf5 in the right format.
All individuals can be sbatched via a script in `packagesSupport/cluster_runs/Antonio_callROH/`
Singe Individuals can be rerun with `Notebooks/ParallelRuns/parallel_antonio19.ipynb`. In this file, there is also the code to combine all indivdiual
outputs into one summary .csv, that can than be used for downstream analysis (such as plotting all individuals)


## Plotting of Overall Results
Key plotting notebooks are found in `notebooks/figures/`
Important Posterior plotting is found in `plot_posetrior.ipynb`.

Geographic and temporal figure production is found in `plot_map_ROH.ipynb`


### Peak Memory Requirements:
For 1x-2x coverage on 1240k, with readcounts from WG (Poisson), such as for much of Antonio2019 data:
Ca. 25gb per Chromosome (during Internal Calculations).
Total runtime single CPU: Ca. 1h.

## Testing BCFTOOLS and PLINK

Notebooks that wrap shell commands for these two are found in `notebooks/PLINK`.

General strategy: Analyze the Mosaic Individuals with copied in ROH blocks. For that run the tools on the datasets with 100 individuals.
First transfer H5 to VCF (for BCFTOOLS also the PL genotype likelihood field in bcftools is needed). Tools for that conversion are found in `PackagesSupport/h5_python/h5_functions.py`.

After running the VCF (saving into output folders, after transforming to rough "HAPSBURG" ROH format), 2) split up the output .csv (or dataframe) into the individual output folders. These can then be analyzed in the same way as HAPSBURG outputs, with tools in `notebooks/Mosaic1000G_Analysis`

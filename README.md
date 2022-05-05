# Development github repository for the software ``hapROH`` and ``hapCON``
Harald Ringbauer, Yilei Huang, May 2022; Code released under GNU General Public License v3.0

This is the code repository intended for contributing developers (i.e. to manage code used in development and publications).  **The release for users is made available as an installable Python software package [hapROH](https://pypi.org/project/hapROH/)**. The official user documentation is available at [readthedocs](https://haproh.readthedocs.io/en/latest/intro.html). It contains quick-start guides as well as usage examples on test data. 


# hapROH
Software to call ROH from ancient and present-day DNA using reference haplotypes.
Author: Harald Ringbauer, September 2020
Code released under GNU General Public License v3.0

Development Code behind python package hapROH. Please refer to our [online documentation](https://haproh.readthedocs.io/en/latest/tutorial.html) for details about the current version and installation.

This is the git repository for development, as well as code used in the [hapROH publication](https://doi.org/10.1038/s41467-021-25289-w).


# hapCon
hapCon is an extension of hapROH for estimaing contamination rate for male aDNA samples.

We have prepared a detailed [online documentation](https://haproh.readthedocs.io/en/latest/hapCON.html) for hapCon. In addition, a jupyter notebook guide for using our method is availble at ./Notebooks/Vignettes/hapCon_vignette.ipynb in this repository. We also have a [preprint](https://www.biorxiv.org/content/10.1101/2021.12.20.473429v1) for hapCon if you are interested in more technical details and its usage limits.

### Install
hapCon is bundled together with hapROH, you only need to install the python package hapROH to use both hapCon and hapROH.

To install,

    pip install hapROH

To upgrade,

    pip install --upgrade hapROH

### Quick Starting Guide
The quickest way to have a test run of hapCon is to use the prepared Python script ./bam/hapCONX.py. It is essentially a wrapper script for the core function of hapCon.

To use the hapCONX.py script, you need at least three input: the pileup file for your sample, the reference panel and the meatadata for the reference panel. You can download the reference panel from https://www.dropbox.com/s/1vv8mz9athedpiq/data.zip?dl=0 (TODO: replace this dropbox link with zenodo later). To generate the pileup file, you can use either [samtools mpileup](http://www.htslib.org/doc/samtools-mpileup.html) or [BamTable](https://bioinf.eva.mpg.de/BamTable/).
    
    python hapCONX.py -m [path to pileup file] -r [path to reference panel] --meta [path to the metadata file]
    
For more details about how to prepare the pileup file and more customized usage of hapCon, please refer to our [jupyter notebook tutorial](https://github.com/hyl317/hapROH/blob/master/Notebooks/Vignettes/hapCon_vignette.ipynb).

###

Author: Yilei Huang, April 2022

# hapROH with Contamination
Joint estimation of ROH blocks and Contamination Rate

Please refer to our [readthedocs](https://haproh.readthedocs.io/en/latest/hapROH_with_contamination.html) site for details.

<!-- # hapROH
Software to call ROH from ancient and present-day DNA using reference haplotypes.
Author: Harald Ringbauer, September 2020
Code released under GNU General Public License v3.0

Development Code behind python package hapROH.
See https://pypi.org/project/hapROH/ for detail about the current version and installation.

This is the git repository for development, as well as code used in the hapROH publication
(Preprint at https://doi.org/10.1101/2020.05.31.126912)

This code is desgined to run on the University of Chicago midway cluster, which uses a `slurm` job scheduling system.

The code behind the hapROH package can be found in `./package/`

Large data is not shared (the `./Data` folder is git-ignored) via github. You will need to download empirical data. You can download the reference data as explained on the Pip Project site. https://pypi.org/project/hapROH/. Alternatively you can download the 1000G data directly and run the notebooks to process them into the required hdf5 format (see below). For simulated data, you can generate data via running the simulation notebooks described below.




# Internal notes about folder structure:

## Data Preparation:
This scripts can prepare the reference data

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

### Produce C extension
In order to run the code, you need to build the C extension that implements the forward/backward pass. If you do not install hapROH via PIP or wish to use the development version, you have to compile the C extension yourself. You can do this via:

Switch to the C folder and build the extensions via:

cd package/hapsburg   
module load gcc  
module load python  
cythonize -a -i cfunc.pyx  
  
This compiles a C file into package/hapsburg. Manually importing `from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind` makes you use these file.d

Tested with gcc/6.1 and python/anaconda-2020.02

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

## \ReichLabEigenstrat
Contains the Eigenstrat from the Reich website. Use script in ./Notebooks/PrepareData/prepare_Reich_Eigenstrat to download and unzip it.



# How to run notebooks on cluster
Assumes that the cluster uses a `slurm` job scheduling system.

First run: 
`sbatch ./jupyter.sbatch`

and then collect the URL from the .err file of the job. 

Afterwards: Do not forget to scancel the job otherwise it runs until the time limit is hit.

This allows one to run memory or cpu-intense tasks via a notebook on the cluster.

For debugging/info: Show job account information for a specific job:
sacct -j jobid --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss,MaxVMSize,nnodes,ncpus,nodelist



## Notebooks to run for creation of data used in paper
Notebooks to run hapROH in parallel on simulated as well as empirical data are found in
`./Notebooks/ParallelRuns/`

these run Hapsburg in a parallelized fashion, assuming that the notebook is run on a cluster with enough nodes (which can be set). The general strategy is to produce lists of parameter files (with one list of parameters per run, e.g. which individual and chromosome), which are then spread out by 'pool.starmap' to nodes via a multirun function.

## Also: hapROH  can be run on cluster with slurm:
Code can be found in ./cluster_runs/

Usually the submission scripts are written as array jobs:
Each folder contains a python file that runs a single individual, taking an integer as input. The scripts map this integer to individuals using the accoring meta files (somewhere in ./data/)

The second file in each folder is `slurm_batch.sh` which is the shell submission script. Modify the .sh scripts to specify what batch to run (i.e. how many and which individuals to submit)

They are submitted with:
`sbatch slurm_batch.sh`

Usually Analysis per chromosomes are run, and results saved into folders: basefolder/iid/chr0/

There is a helper function `combine_individual_data` from `./packagesSupport/parallel_runs/helper_functions` that combines these results per Chromosome into a .csv per Individual, with a option to delete the Chromosome Result folder to save space. For detailled downstream analysis (e.g. of the Posterior), the `delete` argument has to be False to keep the posterior file. The standard output ist `base_folder/IID_roch_full.csv`

## Further Processing per Individual:
These jobs produce summary files for each individual (ususally with `iid_roh_full.csv`) with fields
(Start,End,StartM,EndM,length,lengthM,iid,ch)

These are then further post-process with functions in `./packagesSupport/pp_individual_roh_csvs.py`.

The ususual structure of doing this is:

1) Produce list of paths to .csvs
2) Iterate over these csvs, extract and postprocess (e.g. gap-merging) these
3) These are then merged into one master .csv file with the fields, saved as `combined_roh.csv` with columns
(iid, pop, max_roh, sum_roh, n_roh, lat, lon, age), where some summary details can be merged back in from the master meta dataframe. The function to do so can be found in pp_individual_roh_csvs.py in PackagesSupport. This function is typically called in the ParallelRun notebooks, towards the end

One example of this workflow: `./Notebooks/PrepareData/roh_process.ipynb`

### Combining Different Datasets
These combined csvs are then combined again in 'combine_roh_results.ipynb' into one giant master dataframe. In this function, also "region" columns are defined, which are the basis for downstream processing.

### Example: Antonio 2019 Individuals
There is a notebook in `Notebooks/PrepareData` to prepare the Meta, as well as the hdf5 in the right format.
All individuals can be sbatched via a script in `packagesSupport/cluster_runs/Antonio_callROH/`
Singe Individuals can be rerun with `./Notebooks/ParallelRuns/parallel_antonio19.ipynb`. In this file, there is also the code to combine all indivdiual
outputs into one summary .csv, that can than be used for downstream analysis (such as plotting all individuals)


## Plotting of Overall Results
Notebooks that plot the data can be found in `./Notebooks/figures/`

In particular, posterior plotting is found in `plot_posetrior.ipynb`.

Geographic and temporal figure production is found in `plot_map_ROH.ipynb`


### Peak Memory Requirements:
For pseudo-haploid 1240k data: 2-10 gb per Chromosome (during Internal Calculations), depending on target coverage.

Total runtime single CPU (Intel Xeon E5-2680 v4 (2.40GHz) processor) for one individual (22 Autosomes): Up to 15 minutes.



## Testing BCFTOOLS and PLINK

Notebooks that wrap shell commands for these two are found in `./Notebooks/PLINK`.

General strategy: Analyze the Mosaic Individuals with copied in ROH blocks. For that run the tools on the datasets with 100 individuals.
First transfer H5 to VCF (for BCFTOOLS also the PL genotype likelihood field in bcftools is needed). Tools for that conversion are found in `./PackagesSupport/h5_python/h5_functions.py`.

After running the VCF (saving into output folders, after transforming to rough "HAPSBURG" ROH format), 2) split up the output .csv (or dataframe) into the individual output folders. These can then be analyzed in the same way as HAPSBURG outputs, with tools in `./Notebooks/Mosaic1000G_Analysis` -->

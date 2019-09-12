# HAPSBURG
Software to call ROHs
Author: Harald Ringbauer, 2019
@all rights reserved

python code/notebooks are all in python3

## Installation:
The HMM Python code needs compiling of a C function:

### To create C Extension from .c run command line (within Python3/):
cythonize -a -i cfunc.pyx

Produces some text - what's important: There should be no "compile error" in the output.

On cluster:
load python first before running cythonize:
module load python/3.7.0
(or whatever floats your python3 boat)

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

## Further Processing per Individual:
These jobs produce summary files for each individual (ususally with _roh_full.csv) with fields
(Start,End,StartM,EndM,length,lengthM,iid,ch)

These are then further post-process with functions in `./packagesSupport/pp_individual_roh_csvs.py`.

The ususual structure of doing this is:

1) Produce list of paths to csvs
2) Iterate over these csvs, extract and postprocess (e.g. gap-merging) these (function: )
3) These are then merged into one master .csv file with the fields, saved as `combinedroh.csv` with columns
(iid, pop, max_roh, sum_roh, n_roh, lat, lon, age), where some summary details can be merged back in from the master meta dataframe.

One example of this workflow: `notebooks/PrepareData/roh_process.ipynb`












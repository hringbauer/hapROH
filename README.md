# HAPSBURG
Software to call ROHs
Author: Harald Ringbauer
@all rights reserved

## Installation:
The HMM Python code needs compiling of the C functions:

### To create C Extension from .c run command line (within Python3/):
cythonize -a -i cfunc.pyx

Produces some text - what's important: That there is no "compile error".

On cluster:
load python first before running cythonize:
module load python/3.7.0

## Data Preparation:
For 1240k reference data:
There is a notebook prepare_1000genomes. It downloads 1000 Genomes to local machine and prepares
downsampled and clean 1240k data, and saves it to hdf5 (for the moment)

Goal: 
Have a vcf - and specifiy which indivdiuals to test
And have a reference vcf - and specifiy which individual to use from there 

## Running the HMM

### Running Modes:
cython=
0: Python Fwd/Bwd Algorithm
1: Full Cython Algorithm
2: Optimized Cython Algorithm (linear Nr. references, full Linkage Map Model)


# Folders:
## Data
./Data Folder is git-ignored. That's where the Data goes into.

## /1000 Genome Folder Structure

### Individuals
Contains one column lists of individuals to use

### Markers
Contains 1240k Eigenstrat .snp, that specifies genomic position and map position of every SNP

### Autosome VCF
Contains the input VCFs, downloaded from 1000 Genomes

### Autosome VCF/Subset
Contains also the processed VCFs, downsampled to 

## HDF5
### \FULLHDF5
Contains the full HDF5s.

### \1240kHDF5
Contains the final product. Downsampled hdf5s to individuals as well 

## /ReichLabEigenstrat
Contains the Eigenstrat from the Reich website. Use script in ./Notebooks/PrepareData/prepare_Reich_Eigenstrat to download and unzip it

## Scripts
Contains scripts outside the main Core engine

### \cluster_runs
Contains Python scripts that can be sbatched to cluster. Some of them written as array jobs. The folder also contains the shell scripts for the array jobs





### To preprocess 1000 Genomes Data for X:
`plink --vcf ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz --extract variants1240k --keep-fam EUR_fam.csv --recode vcf --out 1000gX1240kEur --biallelic-only strict --keep-allele-order`

Requires EXACTLY Plink 1.9.

# Diverse

## Run Notebooks on cluster:
sbatch ./jupyter.sbatch

and then collect the URL from the .err file of the job. Afterwards: Do not forget to scancel the job!!


This allows one to run memory or cpu-intense tasks via a notebook on the cluster.

For debugging/info: Show job account information for a specific job:
sacct -j jobid --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss,MaxVMSize,nnodes,ncpus,nodelist









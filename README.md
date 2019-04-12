# HAPSBURG
Software to call ROHs
Author: Harald Ringbauer
@all rights reserved

### To create C Extension from .c run command line:
cythonize -a -i cfunc.pyx

Produces some text - what's important that there is no "compile error".



### Running Modes:
cython=
0: Python Fwd/Bwd Algorithm
1: Full Cython Algorithm
2: Optimized Cython Algorithm (linear Nr. references, full Linkage Map Model)


### Data:
./Data Folder is git-ignored. That's where the Data goes into.

### To preprocess 1000 Genomes Data for X:
`plink --vcf ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz --extract variants1240k --keep-fam EUR_fam.csv --recode vcf --out 1000gX1240kEur --biallelic-only strict --keep-allele-order`

Requires EXACTLY Plink 1.9.


# 1000 Genome Folder Structure

### Individuals
Contains one column lists of individuals to use

### Markers
Contains 1240k Eigenstrat .snp, that specifies genomic position and map position of every SNP

### Autosome VCF
Contains the input VCFs, downloaded from 100 Genomes

### Autosome VCF/Subset
Contains also the processed VCFs, downsampled to 

### HDF5
### \FULLHDF5
Contains the full HDF5s.

### \1240kHDF5
Contains the final product. Downsampled hdf5s to individuals as well 

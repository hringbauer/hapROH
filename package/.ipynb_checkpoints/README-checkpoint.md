# HAPSBURG
Software to call ROHs
Author: Harald Ringbauer, 2020

This package contains functions and wrappers to call ROH from 1240k data, as well as visualizing the results.

### Installation
Please install the package using the Package manager pip. I distribute the source, for building the c extension, see below.

### Example Use
Please find an example notebook, walking through a typical usecase, at XXX.

### c Extension
The package is distributed via source, please build the c extension yourself (ideally via the package cython, set CYTHON=True in setup.py).
The heavy lifting is coded into a cfunction cfunc.c, that was built with cython from cfunc.pyx

### Scope
Standard parameters are tuned for human 1240k capture data (or downsampled SNPs) with 1000 Genome haplotypes as reference, and the software worked for a wide range of test cases. In the first version, HAPSBURG works on eigenstrat file, a future release will add functionality to use diploid genotype calls, or read counts from a .vcf.

If you have whole genome data available, downsample to biallelic 1240k SNPs first.

In case of applications to other kind of SNP or bigger SNP sets, or other organisms, you can contact me at:

harald_ringbauer AT hms harvard edu
(fill in blanks with dots)

Please cite as
[TO BE ANNOUNCED]



# HAPSBURG
Software to call ROHs
Author: Harald Ringbauer, 2020

This package contains functions and wrappers to call ROH from 1240k data using a modern reference panel, and functions to visualize the results.

### Installation
Please install the package using the Package manager pip.

```
python3 -m pip install --user --index-url https://test.pypi.org/simple/ --no-deps hapsburg
```

[IN FINAL VERSION, no --index-url necessary]

I distribute the source only, the setup.py contains information that should help automatically building the necessary c extension.f

For more info about building the c extension, see below. If you install via **pip**, then 

### Scope
Standard parameters are tuned for human 1240k capture data (or downsampled SNPs) with 1000 Genome haplotypes as reference, and the software worked for a wide range of test cases. In the first version, HAPSBURG works on eigenstrat file, a future release will add functionality to use diploid genotype calls, or read counts from a .vcf.

If you have whole genome data available, downsample to biallelic 1240k SNPs first.

In case you are planning applications to other kind of SNP or bigger SNP sets, or other organisms, you can contact me at:

harald_ringbauer AT hms harvard edu
(fill in blanks with dots)


### Get reference Data
Hapsburg currently uses 1000G haplotypes (n=5008), filtered down to bi-allelic 1240k SNPs, including a genetic map. 
We use .hdf5 format for the final output
You can download the prepared reference data (including a necessary metadata .csv) from:  
https://www.dropbox.com/s/0qhjgo1npeih0bw/1000g1240khdf5.tar.gz?dl=0

and unpack into a directory of your choise using 

```
tar -xvf FILE.tar.gz
```

You can then set the link to the folder in the HAPSBURG run parameters. 

You can also download some example Eigenstrats:  
https://www.dropbox.com/s/hjthy138c5t8elv/freilich20.tar.gz?dl=0


### Example Use
Please find an example notebook, walking through a typical usecase, at

[FILL IN FOLDER OF EXAMPLE NOTEBOOK]


### c Extension
The package is distributed via source. This means the c extension has to be built. Ideally, this is done automatically via the package cython, set CYTHON=True in setup.py.

The heavy lifting is coded into a cfunction cfunc.c, that was built with cython from cfunc.pyx

You can also set CYTHON=FALSE, then the extension is compied from cfunc.c directly.


## Citation

Please cite as
[TO BE ANNOUNCED]






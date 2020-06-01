# hapROH
Software to call ROHs in ancient and present-day DNA, using a haplotype reference panel.
Author: Harald Ringbauer, 2020

This package contains functions and wrappers to call ROH for ancient DNA data (for data on 1240k SNPs) and functions to visualize the results.
It needs a haplotype reference panel, 

### Installation
Youc can install the package using the Package manager pip:

```
python3 -m pip install hapROH
```
(`python3 -m` makes sure you use your python installation)


The package distributes source code. The setup.py contains information that should automatically build the necessary c extension.
In case you want to manually build the c extension, find more info in the section below (`c Extension`)

### Scope of the Method
Standard parameters are tuned for human 1240k capture data (1.2 million SNPs) and using 1000 Genome haplotypes as reference. The software worked for a wide range of test cases, both 1240k data and also whole genome sequencing data downsampled to 1240k. Test cases included 45k year old Ust Ishim man, and both American, Eurasian and Oceanian ancient DNA, showing that the method generally works for split times of reference panel and target up to a few 10k years (Neanderthals and Denisovans do not fall into that range).

In the first version, hapROH works on eigenstrat file (either packed or unpacked, the mode can be set). A future release will add functionality to use diploid genotype calls, or genotype likelihoods from a .vcf.

If you have whole genome data available, you should downsample to biallelic 1240k SNPs first.

In case you are planning applications to other kind of SNP or bigger SNP sets, or other organisms, the method parameters have to be updated. You can mirror our procedure (described in the publication), or contact me for assistance at:

harald_ringbauer AT hms harvard edu
(fill in blanks with dots)


### Get reference Data
Hapsburg currently uses global 1000 Genome haplotypes (n=5008), filtered down to bi-allelic 1240k SNPs, including a genetic map. 
We use .hdf5 format for the final output
You can download the prepared reference data (including a necessary metadata .csv) from:  
https://www.dropbox.com/s/0qhjgo1npeih0bw/1000g1240khdf5.tar.gz?dl=0

and unpack into a directory of your choise using 

```
tar -xvf FILE.tar.gz
```

You can then set the link to the folder in the hapROH run parameters. 
You can also download some example Eigenstrats:  
[WILL BE FILLED IN]


### Example Use
Please find an example notebook, walking through a typical usecase, at

./Notebooks/test_pypi_package.ipynb [TEMPORARY, WILL BE FILLED IN]

### Package Requirements
The basic requirements for calling ROH are rather minimal. There are extra Python packages that you need to install (e.g. via `pip`) if you want to use the plotting functionality, you should have `matplotlib` installed. For plotting of maps, you will need `basemap` (warning: installing can be tricky on some systems). If you want to use the effective population size fitting functionality from ROH output, you will need to have `statsmodels` installed.


### c Extension
The heavy lifting of the algorithm is coded into a cfunction cfunc.c, which is built via cython from cfunc.pyx

The package is distributed via source. This means the c extension has to be built. Ideally, this is done automatically via the package cython, as CYTHON=True in setup.py by default.

You can also set CYTHON=FALSE, then the extension is compiled from cfunc.c directly (not tested on all platforms)


## Citation
If you want to cite this software:
[TO BE ANNOUNCED]







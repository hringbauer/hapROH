# hapROH
A program to identify runs of homozygosity (ROH) in ancient and present-day DNA, using a panel of reference haplotypes. It works on data for 1240K SNPs and in the `eigenstrat` format (widely used in ancient DNA). Recommended is at least 0.3x coverage - and the software works assuming pseudo-haploid data.

This package contains functions and wrappers to call ROH and functions for downstream analysis of the results (visualization and analysis). 

For downward compatibility, the package uses `hapsburg` as module name, after installation you can import functions via
`from hapsburg.XX import YY`

## Installation
You can install the package using the Package manager pip:

```
python3 -m pip install hapROH
```
(`python3 -m` makes sure you use your python installation)

If you have it already installed and want to upgrade to a newer hapROH version you can use:
```
python3 -m pip install hapROH --upgrade
```

The package distributes source code. The setup.py contains information that should automatically install the package.
For customized installations, find more info in the section below (`c Extension`)

## Getting Started
To get started, please find vignette jupyter notebooks:
https://www.dropbox.com/sh/eq4drs62tu6wuob/AABM41qAErmI2S3iypAV-j2da?dl=0

These are a ressource to do show example usecases, that you can use as template for your own applications.

These notebooks walk you through examples for 
1) how to use the core functions to call ROH from eigenstrat files, and generate ROH tables from results of multiple individuals ('callROH_vignette')
2) how to use functions for visualizing ROH results ('plotting_vignette' - warning: Some of these are experimental and require additional packages. You might want to consider creating your own plotting functions for visualizing the results in the way that works best for you)
3) how to call IBD on the X chromosome between two male X chromosomes ('callIBD_maleX_vignette', warning: experimental)

## Scope of the Method
Standard parameters are tuned for human 1240K capture data (ca. 1.2 million SNPs used widely in human aDNA analysis) and using 1000 Genome haplotypes as reference. The software is tested on a wide range of test applications, both 1240K data and also whole genome sequencing data downsampled to 1240K SNPs. Successful cases include 45k year old Ust Ishim man, and a wide range of American, Eurasian and Oceanian ancient DNA, showing that the method generally works for split times of reference panel and target up to a few 10k years, which includes all out-of-Africa populations (Attention: Neanderthals and Denisovans do not fall into that range, additionally some Subsaharan hunter gatherer test cases did not give satisfactory results).

In this version, hapROH works on eigenstrat file (either packed or unpacked, the mode can be set). A planned future release will add functionality to use genotype likelihoods from a .vcf.

If you have whole genome data available, you should downsample an create eigenstrat files for biallelic 1240k SNPs first.

In case you are planning applications to other kind of SNP or bigger SNP sets, or even other organisms, the method parameters have to be adjusted (the default parameters are specifically optimized for human 1240K data). You can mirror our procedure to find good parameters (described in the publication), and if you contact me for assistance - I am happy to share my own experience.


## Download reference Data
hapROH currently uses global 1000 Genome data (n=5008 haplotypes), filtered down to bi-allelic 1240K SNPs. 
We use .hdf5 format for the reference panel - which includes a genetic map.

You can download the prepared reference data (including a necessary metadata .csv) from:  
https://www.dropbox.com/s/0qhjgo1npeih0bw/1000g1240khdf5.tar.gz?dl=0

and unpack it using 

```
tar -xvf FILE.tar.gz
```

You then have to link the paths in the hapROH run parameters (see vignette notebook)


## Example Use Case: Vignettes
Please find example notebooks, walking you through a typical application to an eigenstrat file at
https://www.dropbox.com/sh/eq4drs62tu6wuob/AABM41qAErmI2S3iypAV-j2da?dl=0

All you need is a Eigenstrat file, and the reference genome data (see link above), and you are good to go to run your own ROH calling!

There is a vignette notebook for...
1) walking you through the calling of ROH (callROH)
2) producing various figures from the output (plotROH)
3) describing the experimental functionality to identify IBD segements between pairs of male X chromosomes (callIBD_maleX)
4) estimating population sizes from inferred ROH, using a likelihood framework (estimateNe)


## Dependencies
The basic requirements for calling ROH are kept minimal and only sufficient for the core ROH calling ('numpy', 'pandas', 'scipy' & 'h5py'). If you want to use extended analysis and plotting functionality: There are extra Python packages that you need to install (e.g. via `pip` or `conda`). 

1) If you want to use the advanced plotting functionality, you need `matplotlib` installed.
2) For plotting of maps, you will need `basemap` (warning: installing can be tricky on some architectures). 
3) If you want to use the effective population size fitting functionality from ROH output, you require the package `statsmodels`.

## c Extension
For performance reasons, the heavy lifting of the algorithm is coded into a c method (cfunc.c). This "extension" is built via cython from cfunc.pyx This should be done automatically via the package cython (as CYTHON=True in setup.py by default).

You can also set CYTHON=False, then the extension is compiled from cfunc.c directly (experimental, not tested on all platforms).

## Development
The code used to develop this package is deposited at the github repository: https://github.com/hringbauer/hapROH

The package is packed in the folder `./package/`. In addition, there are a large number of notebooks used to test and extensively use the functionality in `./notebooks/`.

## Citation
If you use the software for a scientific publication and want to cite it, you can use:
https://www.biorxiv.org/content/10.1101/2020.05.31.126912v1

## Contact
If you have bug reports, suggeestions or general comments, please feel always free to contact me. I am happy to hear from you. Bug reports and user suggestions will help me to improve this software - so please do not hesitate to reach out!

harald_ringbauer AT hms harvard edu
(fill in blanks with dots)

## 	Acknowledgments
Big thank you to my co-authors Matthias Steinrücken and John Novembre. The project profited immensely from Matthias' deep knowledge about HMMs and from John's extensive experience in developing population genetics software. Countless discussions with both were key for moving forward this project. Another big thanks goes to Nick Patterson, who informed me about the benefits of working with rescaled HMMs - substantially improving the runtime of hapROH. I also want to thank users who find and report software bugs (Mélanie Pruvost, Ke Wang). Their feedback helped to identify and remove errors in the program.

Author:
Harald Ringbauer, 2020
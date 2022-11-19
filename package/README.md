# hapROH & hapCon
This Python package contains two softwares:

1) hapROH
This software identifies runs of homozygosity (ROH) in ancient and present-day DNA by using a panel of reference haplotypes. This package contains functions and wrappers to call ROH and functions for downstream analysis and visualization.

For downward compatibility, the package uses `hapsburg` as module name. After installation you can import Python functions via
`from hapsburg.XX import YY`

2) hapCon
This software estimates contamination in male X chromosome via using a panel of reference haplotypes. It has been incorporate into the hapROH package since version 0.4a1, no additional installation is needed. It has the same dependencies as hapROH. 

Detailled Instructions for Installation, Getting Started and Vignettes can be found at:
https://haproh.readthedocs.io


## Scope of the Method

### hapROH
Standard parameters are tuned for human 1240K capture data (ca. 1.2 million SNPs used widely in human aDNA analysis) and using 1000 Genome haplotypes as reference. The software is tested on a wide range of test applications, both 1240K data and also whole genome sequencing data downsampled to 1240K SNPs. Successful cases include 45k year old Ust Ishim man, and a wide range of American, Eurasian and Oceanian ancient DNA, showing that the method generally works for split times of reference panel and target up to a few 10k years, which includes all out-of-Africa populations (Attention: Neanderthals and Denisovans do not fall into that range, additionally some Subsaharan hunter gatherer test cases did not give satisfactory results).

Currently, hapROH works on data for 1240K SNPs and in unpacked or packed `eigenstrat` format (which is widely used in human ancient DNA). The software assumes pseudo-haploid or diploid genotype data (the mode can be set, by default it is pseudo-haploid). The recommended coverage range is 400,000 or more 1240K SNPs covered at least once (i.e. at least ca. 0.3x coverage).

If you have whole genome data available, you have to downsample an create eigenstrat files for biallelic 1240k SNPs first.

In case you are planning applications to other kind of SNP or bigger SNP sets, or even other organisms, the method parameters have to be adjusted (the default parameters are specifically optimized for human 1240K data). You can mirror our procedure to find good parameters (described in the publication), and if you contact me for assistance - I am happy to share my own experience.

### hapCon
It works directly from BAM file or from samtools mpileup output. We have created two reference panels for hapCON: one for 1240k data and the other for WGS data. Standard parameters are tuned towards these two use cases.


## Dependencies
The basic requirements for calling ROH are kept minimal and only sufficient for the core ROH calling ('numpy', 'pandas', 'scipy', 'numdifftools', 'h5py'). If you want to use extended analysis and plotting functionality: There are extra Python packages that you need to install (e.g. via `pip` or `conda`). 

1) If you want to use the advanced plotting functionality, you need `matplotlib` installed.
2) For plotting of maps, you will need `basemap` (warning: installing can be tricky on some architectures as C packages are required). 
3) If you want to use the effective population size fitting functionality from ROH output, you require the package `statsmodels`.

## c Extension
For performance reasons, the heavy lifting of the algorithm is coded into a c method (cfunc.c). This "extension" is built via cython from cfunc.pyx This should be done automatically via the package cython (as CYTHON=True in setup.py by default).

You can also set CYTHON=False, then the extension is compiled from cfunc.c directly (experimental, not tested on all platforms).

## Development
The code used to develop this package is deposited at the github repository: 
https://github.com/hringbauer/hapROH

The package is packed in the folder `./package/`. In addition, there are a large number of notebooks used to test and extensively use the functionality in `./notebooks/`.

## Citation
If you use the software of hapROH for a scientific publication and want to cite it, you can use:
https://www.biorxiv.org/content/10.1101/2020.05.31.126912v1

## Contact
If you have bug reports, suggestions or general comments, please contact me. I am happy to hear from you. Bug reports and user suggestions will help me to improve this software - so please do not hesitate to reach out!

harald_ringbauer AT eva mpg de
(fill in blanks with dots and AT with @)

## 	Acknowledgments
Big thank you to the original co-authors Matthias Steinrücken and John Novembre. The project profited immensely from Matthias' deep knowledge about HMMs and from John's extensive experience in developing population genetics software. Countless discussions with both have been key for moving forward this project. Another big thanks goes to Nick Patterson, who informed me about the benefits of working with rescaled HMMs - substantially improving the runtime of hapROH. 

I want to acknowledge users who find and report software bugs (Mélanie Pruvost, Ke Wang, Ruoyun Hui, Selina Carlhoff, Matthew Mah, Xiaowen Jia) and all users who reached out with general questions and requests (Rosa Fregel, Federico Sanchez). This feedback has helped to remove errors in the program and to improve its usability. Many thanks!



Authors:
Harald Ringbauer, Yilei Huang, 2021
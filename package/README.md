# hapROH & hapCon
This Python package contains two computational tools for ancient DNA, ***hapROH*** and ***hapCon***.

**For detailed Instructions for Installation and Manuals, including Getting Started and Vignettes, please visit the official tutorial:
https://haproh.readthedocs.io**

1) hapROH
This software identifies runs of homozygosity (ROH) in ancient and present-day DNA by using a panel of reference haplotypes. This package contains functions and wrappers for calling ROH, as well as functions for downstream analysis and visualization.

For downward compatibility, the package uses `hapsburg` as the module name. After installation, you can import Python functions via
`from hapsburg.XX import YY`

2) hapCon
This software estimates contamination in the male X chromosome by using a panel of reference haplotypes. It has been incorporated into the hapROH package since version 0.4a1; no additional installation is needed. 


## Scope

### hapROH
Standard parameters are tuned for human 1240K capture data (approximately 1.2 million SNPs widely used in human aDNA analysis) and utilize 1000 Genomes haplotypes as a reference. The software is tested on a wide range of test applications, including both 1240K data and whole-genome sequencing data downsampled to 1240K SNPs. Successful cases include 45k year old Ust Ishim man, and a wide range of American, Eurasian and Oceanian ancient DNA, showing that the method generally works for split times of reference panel and target up to a few 10k years, which includes all out-of-Africa populations (Attention: Neanderthals and Denisovans do not fall into that range, additionally some Subsaharan hunter gatherer test cases did not give satisfactory results).

Currently, hapROH works on data for 1240K SNPs and in unpacked or packed `eigenstrat` format (which is widely used in human ancient DNA). The software assumes pseudo-haploid or diploid genotype data (the mode can be set; by default, it is pseudo-haploid). The recommended coverage range is 400,000 or more 1240K SNPs covered at least once (i.e., at least >0.3x coverage).

If you have whole-genome data available, you must first downsample and create eigenstrat files for the biallelic 1240k SNPs.

In case you plan to apply this method to other types of SNPs, larger SNP sets, or other organisms, the method parameters must be adjusted (the default parameters are explicitly optimized for human 1240K data). You can mirror our procedure to find good parameters (described in the publication), and if you contact me for assistance, I am happy to share my own experience.

### hapCon
This software works directly from a BAM file or from samtools mpileup output. We have created two reference panels for hapCON: one for 1240k data and the other for WGS data. The standard parameters are tuned towards these two use cases.

## Updates:
The text file`./change_log.md` describes updates in the various versions of this software.


## Dependencies
The basic requirements for calling ROH are kept minimal, allowing only the core ROH calling ('numpy', 'pandas', 'scipy', 'numdifftools', 'h5py'). If you want to use extended analysis and plotting functionality. There are extra Python packages that you need to install (e.g., via `pip` or `conda`). 

1) If you want to use the advanced plotting functionality, you need `matplotlib` 
2) For plotting maps, you will need `basemap` (warning: installing can be tricky on some architectures as C packages are required). 
3) If you want to use the effective population size fitting functionality from ROH output, you require the package `statsmodels`.


## c Extension
For performance reasons, the heavy lifting of the algorithm is coded into a c method (cfunc.c). This "extension" is built via cython from `cfunc.pyx` This should be done automatically via the package cython (as CYTHON=True in setup.py by default).

You can also set `CYTHON=False`, then the extension is compiled from `cfunc.c` directly (experimental, not tested on all platforms).


## Software Development
The code used to develop this package is deposited at the github repository: 
https://github.com/hringbauer/hapROH

The package is packed in the folder `./package/`. In addition, there are a large number of notebooks used to test and extensively use the functionality in `./notebooks/`.


## Citation
If you use the software for a scientific publication and want to cite it, you can use:
- **hapROH***: https://doi.org/10.1038/s41467-021-25289-w
- **hapCon***: https://doi.org/10.1093/bioinformatics/btac390


## Contact
If you have bug reports, suggestions or general comments, please contact us. We are happy to hear from you. Bug reports and user suggestions will help me to improve this software - so please do not hesitate to reach out!

harald_ringbauer AT eva mpg de
yilei_huang AT eva.mpg.de

(fill in blanks with dots and AT with @)


## 	Acknowledgments
Big thank you to the original co-authors Matthias Steinrücken and John Novembre. The project profited immensely from Matthias' deep knowledge about HMMs and from John's extensive experience in developing population genetics software. Countless discussions with both have been key for moving forward this project. Another big thanks goes to Nick Patterson, who informed me about the benefits of working with rescaled HMMs - substantially improving the runtime of hapROH. 

I want to acknowledge users who find and report software bugs (Mélanie Pruvost, Ke Wang, Ruoyun Hui, Selina Carlhoff, Matthew Mah, Xiaowen Jia) and all users who reached out with general questions and requests (Rosa Fregel, Federico Sanchez). This feedback has helped to remove errors in the program and to improve its usability. Many thanks!


Authors:
Harald Ringbauer, Yilei Huang, 2024
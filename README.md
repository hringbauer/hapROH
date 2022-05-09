# Development github repository for the software ``hapROH`` and ``hapCON``
Harald Ringbauer, Yilei Huang, May 2022; Code released under GNU General Public License v3.0

**The release for users is made available as an installable Python software package [hapROH](https://pypi.org/project/hapROH/)**. The **official user documentation is available at [readthedocs](https://haproh.readthedocs.io/en/latest/intro.html)**. It contains quick-start guides as well as usage examples on test data. This here is the github code repository intended for contributing developers (i.e. to manage code used in development and publications). 

Scientific manuscripts describing the methods are available at:
- https://doi.org/10.1038/s41467-021-25289-w for hapROH
- https://doi.org/10.1101/2021.12.20.473429 for hapCon



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

# Development github repository for the software ``hapROH`` and ``hapCON``
Harald Ringbauer, Yilei Huang, May 2022; Code released under GNU General Public License v3.0

This is the code repository intended for contributing developers (i.e. to manage code used in development and publications).  **The release for users is made available as an installable Python software package [hapROH](https://pypi.org/project/hapROH/)**. The official user documentation is available at [readthedocs](https://haproh.readthedocs.io/en/latest/intro.html). It contains quick-start guides as well as usage examples on test data. 

There are three main parts of the hapROH package:

## hapROH
Software to call ROH from ancient and present-day DNA using reference haplotypes.
Author: Harald Ringbauer, September 2020
Code released under GNU General Public License v3.0

Development Code behind the Python package hapROH. Please refer to our [online documentation](https://haproh.readthedocs.io/en/latest/tutorial.html) for details about the current version and installation.

This is the git repository for development, as well as code used in the [hapROH publication](https://doi.org/10.1038/s41467-021-25289-w).


## hapCon
hapCon is an extension of hapROH for estimating contamination rate for male aDNA samples.

We have prepared a detailed [online documentation](https://haproh.readthedocs.io/en/latest/hapCON.html) for hapCon. In addition, a Jupyter notebook guide for using our method is available at ./Notebooks/Vignettes/hapCon_vignette.ipynb in this repository. We also have a [preprint](https://www.biorxiv.org/content/10.1101/2021.12.20.473429v1) for hapCon if you are interested in more technical details and its usage limits.

## hapROH with Contamination
Joint estimation of ROH blocks and Contamination Rate

Please refer to our [readthedocs](https://haproh.readthedocs.io/en/latest/hapROH_with_contamination.html) site for details.

# Installation
You only need to install the PyPi package hapROH to use both hapCon and hapROH:

To install via pip:

    pip install hapROH

To upgrade an existing installation:

    pip install --upgrade hapROH

# References
Scientific manuscripts introducing the methods are available at:
- https://doi.org/10.1038/s41467-021-25289-w for hapROH
- https://doi.org/10.1101/2021.12.20.473429 for hapCon

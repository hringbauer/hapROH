hapCON
==========================================================================


This software can estimate contamination in male X chromosome using a haplotype copying framework. Details of the model are described in a `preprint <https://doi.org/10.1101/2021.12.20.473429 >`_ 

An implementation (hapCON) has been incorporated into the hapROH package since version 0.4a1. No additional installation is needed. hapCON works directly from BAM file or from samtools mpileup output. We have created two reference panels for common use cases in human aDNA data: One for 1240k data and the other for WGS data.

The core functionality of hapCON is exposed via :meth:`hapsburg.PackagesSupport.hapsburg_run.hapCon_chrom_BFGS`.

A short tutorial and example usage is available `here <https://github.com/hyl317/hapROH/blob/master/Notebooks/Vignettes/hapCON_vignette.ipynb>`_.

Authors: Yilei Huang, Harald Ringbauer Dec 2021

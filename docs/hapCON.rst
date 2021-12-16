hapCON: An Extension of hapROH for Estimating Contamination in maleX
==========================================================================


hapCON is a new method for estimating contamination in male X chromosome. It has been incorporate into hapROH since version xxx, so no additional installation is needed. It has the same dependencies as hapROH. It works directly from BAM file or from samtools mpileup output. We have created two reference panels for hapCON: one for 1240k data and the other for WGS data.

The core functionality of hapCON is exposed via :meth:`hapsburg.PackagesSupport.hapsburg_run.hapCon_chrom_BFGS`.

A short tutorial and example usage can be found `here <https://github.com/hyl317/hapROH/blob/master/Notebooks/Vignettes/hapCON_vignette.ipynb>`_.

Author: Yilei Huang, Dec 2021
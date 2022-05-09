hapCon
==========================================================================

Scope of the Method
**************************
This software can estimate contamination in male X chromosome using a haplotype copying framework. The model is described in a `preprint <https://doi.org/10.1101/2021.12.20.473429>`_.

Application Range of hapCon:
    HapCon works with ancient DNA data for modern human. It does not work for Neaderthals and Denisovans and other archaic hominins. 
    Our simulations show that hapCon works well for a wide range of ancestries and sample age, including Ust Ishim, one of the oldest modern humans ever sequenced, ancient Africans (e.g, Mota) and native Americans.
    However, we have found that hapCon does not work well with Sub-Saharan ancient foragers (e.g, samples from Site Hora 1 and Fingira, both in present-day Malawi, `Lipson et al. Nature 2022 <https://www.nature.com/articles/s41586-022-04430-9>`_). These samples can contain a substantial amount of south African related ancestry (such as represented by present-day groups—Juǀ'hoansi (San)), which is not present in the 1000Genome dataset. For details, see Supplementary Note 4 of our manuscript.


Getting started
*************************

hapCon has been incorporated into the hapROH package (since version 0.4a1). For example use cases, that you can adapt for your purpose, please see our `vignette notebook <https://github.com/hyl317/hapROH/blob/master/Notebooks/Vignettes/hapCon_vignette.ipynb>`_.


The quickest way to have a test run of hapCon is to use the prepared Python script ./bam/hapCONX.py. It is a wrapper script for the core function of hapCon.

To use the hapCONX.py script, you need at least three input: the pileup file for your sample, the reference panel and the metadata for the reference panel. You can download the reference panel from https://www.dropbox.com/s/1vv8mz9athedpiq/data.zip?dl=0 (TODO: replace this dropbox link with zenodo later). To generate the pileup file, you can use either `samtools mpileup <http://www.htslib.org/doc/samtools-mpileup.html>`_ or `BamTable <https://bioinf.eva.mpg.de/BamTable/>`_.
    
.. code-block:: console

    python hapCONX.py -m [path to pileup file] -r [path to reference panel] --meta [path to the metadata file]
    
    
For more details about how to prepare the pileup file and more customized usage of hapCon, please check out the `jupyter notebook tutorial <https://github.com/hyl317/hapROH/blob/master/Notebooks/Vignettes/hapCon_vignette.ipynb>`_.


Input
*************************

hapCon works directly from BAM file or from `samtools mpileup <http://www.htslib.org/doc/samtools-mpileup.html>`_ or `BamTable <https://bioinf.eva.mpg.de/BamTable/>`_ output. 

We have created two reference panels for common use cases in human aDNA data: One for 1240k data and the other for WGS data (TODO: add a link to zenodo repo after paper acceptance).

The core functionality of hapCon is exposed via :meth:`hapsburg.PackagesSupport.hapsburg_run.hapCon_chrom_BFGS`. The input can be  
1) BAM file
2) output from `samtools mpileup <http://www.htslib.org/doc/samtools-mpileup.html>`_ 
3) output from `BamTable <https://bioinf.eva.mpg.de/BamTable/>`_. 

We recommend using BamTable for preprocessing your BAM file as it provides the most flexibility. 

Input BAM File
*************************
The input BAM file should have passed common preprocessing steps depending on your aDNA data type (e.g. which UDG treatment, double or single stranded library preparation protocol, 1240k capture or shotgun data). It is the same preprocessing steps as done for genotype calling and includes merging paried reads, PCR deduplication, filtering to read lengths and qualities. For a set of commonly used preprocessing steps in aDNA research we refer to `Eager 2  <https://github.com/nf-core/eager>`.
    
Base qaulity or alignment quality filtering can be done during the process of generating read counts with `samtools <http://www.htslib.org/doc/samtools.html>`_ or `BamTable <https://bioinf.eva.mpg.de/BamTable/>`_.


Usage Notes:
    
    * If your data contains African ancestry, please adjust the parameter "exclude_pops" in :meth:`hapsburg.PackagesSupport.hapsburg_run.hapCon_chrom_BFGS`. As explained in our manuscript, by default it excludes African haplotypes in the reference panel as this alleviates the "attraction effect". In case when the sample has African ancestry, however, the whole reference panel should be used.


Example Use Case: Vignettes
*****************************
For example use cases, please checkout our `tutorial <https://github.com/hyl317/hapROH/blob/master/Notebooks/Vignettes/hapCon_vignette.ipynb>`_.


Authors: Yilei Huang, Harald Ringbauer April 2022

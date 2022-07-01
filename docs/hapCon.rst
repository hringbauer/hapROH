hapCon
==========================================================================

This software estimates contamination of human aDNA data using a haplotype copying framework for male X chromosomes. The tool is published in `Bioinformatics <https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac390/6607584>`_.


Application range of hapCon
**************************
HapCon works with ancient DNA data for modern human. It does not work for Neaderthals and Denisovans and other archaic hominins.  Our simulations show that hapCon works well for a wide range of ancestries and sample age, including Ust Ishim, one of the oldest modern humans ever sequenced, ancient Africans (e.g, Mota) and native Americans. However, we have found that hapCon does not work well with Sub-Saharan ancient foragers (e.g, samples from Site Hora 1 and Fingira, both in present-day Malawi, `Lipson et al. Nature 2022 <https://www.nature.com/articles/s41586-022-04430-9>`_). These samples can contain a substantial amount of south African related ancestry (such as represented by present-day groups—Juǀ'hoansi (San)), which is not present in the 1000Genome dataset. For details, see Supplementary Note 4 of our manuscript.


Input
*************************

hapCon works directly from BAM file or from `samtools mpileup <http://www.htslib.org/doc/samtools-mpileup.html>`_ or `BamTable <https://bioinf.eva.mpg.de/BamTable/>`_ output. 

We have created two reference panels for common use cases in human aDNA data: One for 1240k data and the other for WGS data (https://doi.org/10.5281/zenodo.6619138).

The core functionality of hapCon is exposed via :meth:`hapsburg.PackagesSupport.hapsburg_run.hapCon_chrom_BFGS`. The input can be any one of the following,

* BAM file
* output from `samtools mpileup <http://www.htslib.org/doc/samtools-mpileup.html>`_ 
* output from `BamTable <https://bioinf.eva.mpg.de/BamTable/>`_. 

We recommend using BamTable for preprocessing your BAM file as it provides the most flexibility. For optimal performance, the input BAM file should have passed common preprocessing steps as required for producing genotype data for your data type (e.g. depending on UDG treatment, double or single stranded library preparation protocol, 1240k capture or shotgun data). This includes merging paired reads, PCR deduplication, filtering by read lengths and mapping qualities. For a set of commonly used QC steps we refer to `Eager2 <https://github.com/nf-core/eager>`_. Base qaulity or alignment quality filtering can be done within `samtools <http://www.htslib.org/doc/samtools.html>`_ or `BamTable <https://bioinf.eva.mpg.de/BamTable/>`_.


Getting started
*************************

hapCon has been incorporated into the hapROH package (since version 0.4a1). For example use cases, that you can adapt for your purpose, please see our `vignette notebook <https://github.com/hyl317/hapROH/blob/master/Notebooks/Vignettes/hapCon_vignette.ipynb>`_.


The quickest way to have a test run of hapCon is to use the prepared Python script ./bam/hapCONX.py. It is a wrapper script for the core function of hapCon.

To use the hapCONX.py script, you need at least three input: the pileup file for your sample, the reference panel and the metadata for the reference panel. 

You can download the reference panel and test BAM file from https://doi.org/10.5281/zenodo.6619138.

Assuming you have the above-mentioned file in your current directory, to generate the pileup file, you can use either `samtools mpileup <http://www.htslib.org/doc/samtools-mpileup.html>`_ or `BamTable <https://bioinf.eva.mpg.de/BamTable/>`_. We will use samtools as an example here,

.. code-block:: console

    samtools mpileup -o test.mpileup -q 30 -Q 30 --positions chrX_1240k.bed SUA001.bam

this generates a pileup file for your BAM. Then we are ready to run hapCon.
    
.. code-block:: console

    python hapCONX.py -m test.mpileup -r chrX_1240k.hdf5 --meta meta_df_all.csv

By the default, the output file will have the same prefix as the .mpileup file, and it will be named as $prefix.hapCon.txt. By using the above input and command, the output will have the following format,

.. code-block:: console

    Number of target sites covered by at least one read: 3999
    Method1: Fixing genotyping error rate
	    Estimated genotyping error via flanking region: 0.001502
	    MLE for contamination using BFGS: 0.102113 (0.076802 - 0.127424)

It gives

#. Number of sites covered by one or more reads (aka the length of the HMM chain)
#. Genotype error estimated by measuring heterozygosity on supposedly non-segregating sites
#. MLE of contamination rate and its confidence intervals

If your data contains African ancestry, please adjust the parameter "--exclude". As explained in our manuscript, by default it excludes African haplotypes in the reference panel as this alleviates the "attraction effect". In case when the sample has African ancestry, however, the whole reference panel should be used. 

If you believe the contamination source of your sample is not of continental European, please adjust --con parameter accordingly. To see a full list of adjustable parameters, 

.. code-block:: console

    python hapCONX.py -h


For more details about the usage of hapCon, please check the next section.


Example Use Case: Vignettes
*****************************
For detailed example use cases, please checkout our `tutorial <https://github.com/hyl317/hapROH/blob/master/Notebooks/Vignettes/hapCon_vignette.ipynb>`_.



Authors: Yilei Huang, Harald Ringbauer May 2022

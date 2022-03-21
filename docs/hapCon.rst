hapCon
==========================================================================


This software can estimate contamination in male X chromosome using a haplotype copying framework. Details of the model are described in a `preprint <https://doi.org/10.1101/2021.12.20.473429>`_.

An implementation (hapCon) has been incorporated into the hapROH package since version 0.4a1. No additional installation is needed.
hapCon works directly from BAM file or from `samtools mpileup <http://www.htslib.org/doc/samtools-mpileup.html>`_ or `BamTable <https://bioinf.eva.mpg.de/BamTable/>`_ output. 
We have created two reference panels for common use cases in human aDNA data: One for 1240k data and the other for WGS data (TODO: add a link to zenodo repo after paper acceptance).

The core functionality of hapCon is exposed via :meth:`hapsburg.PackagesSupport.hapsburg_run.hapCon_chrom_BFGS`. It can take as input BAM file, output from `samtools mpileup <http://www.htslib.org/doc/samtools-mpileup.html>`_ or output from `BamTable <https://bioinf.eva.mpg.de/BamTable/>`_. We recommend using BamTable for preprocessing your BAM file as it provides the most flexibility. For more details, please checkout our `tutorial <https://github.com/hyl317/hapROH/blob/master/Notebooks/Vignettes/hapCon_vignette.ipynb>`_.

Application Range of hapCon:
    HapCon works with ancient DNA data for modern human. It does not work for Neaderthals and Denisovans and other archaic hominins. 
    Our simulations show that hapCon works well for a wide range of ancestries and sample age, including Ust Ishim, one of the oldest modern humans ever sequenced, ancient Africans (e.g, Mota) and native Americans.
    However, we have found that hapCon does not work well with Sub-Saharan ancient foragers (e.g, samples from Site Hora 1 and Fingira, both in present-day Malawi, `Lipson et al. Nature 2022 <https://www.nature.com/articles/s41586-022-04430-9>`_). These samples can contain a substantial amount of south African related ancestry (such as represented by present-day groups—Juǀ'hoansi (San)), which is not present in the 1000Genome dataset. For details, see Supplementary Note 4 of our manuscript.

Usage Notes:
    * For the input BAM file, it should have passed a few common QC steps, like PCR deduplication, filter for read lengths, etc. Base qaulity or alignment quality filtering can be done within samtools or BamTable. We suggest performing the same QC steps (with the exception of PMDtools) for hapCon as you would do for other common popgen analysis, like PCA or F-statistics.
    * If your data contains African ancestry, please adjust the parameter "exclude_pops" in :meth:`hapsburg.PackagesSupport.hapsburg_run.hapCon_chrom_BFGS`. As explained in our manuscript, by default it excludes African haplotypes in the reference panel as this alleviates the "attraction effect". In case when the sample has African ancestry, however, the whole reference panel should be used.


Authors: Yilei Huang, Harald Ringbauer March 2022

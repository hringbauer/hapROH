hapCon
==========================================================================


This software can estimate contamination in male X chromosome using a haplotype copying framework. Details of the model are described in a `preprint <https://doi.org/10.1101/2021.12.20.473429>`_ 

An implementation (hapCon) has been incorporated into the hapROH package since version 0.4a1. No additional installation is needed. hapCon works directly from BAM file or from samtools mpileup output. We have created two reference panels for common use cases in human aDNA data: One for 1240k data and the other for WGS data.

The core functionality of hapCon is exposed via :meth:`hapsburg.PackagesSupport.hapsburg_run.hapCon_chrom_BFGS`.

A short tutorial and example usage is available `here <https://github.com/hyl317/hapROH/blob/master/Notebooks/Vignettes/hapCon_vignette.ipynb>`_.

Application Range of hapCon:
    HapCon works with ancient DNA data for modern human. It does not work for Neaderthals and Denisovans and other archaic hominins. 
    Our simulations show that hapCon works well for a wide range of ancestries and sample age, including Ust Ishim, one of the oldest modern humans ever sequenced, ancient Africans (e.g, Mota) and native Americans.
    However, we have found that hapCon does work well with Sub-Saharan ancient foragers (e.g, samples from Site Hora 1 and Fingira, both in present-day Malawi, `Lipson et al. Nature 2022 <https://www.nature.com/articles/s41586-022-04430-9>`_). These samples can contain a substantial amount of Southern African related ancestry, which is not present in the 1000Genome dataset. For details, see Supplementary Note 4 of our manuscript.


Authors: Yilei Huang, Harald Ringbauer March 2022

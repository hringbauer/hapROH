Frequently Asked Questions (FAQ)
==================================

I have a .bam, what should I do?
**********************************

For ancient DNA data, *hapROH* takes `pseudo-haploid eigenstrat files <https://reich.hms.harvard.edu/software/InputFileFormats>`_ as input. 
As of 2021, this data type is widely used in human ancient DNA data analysis (such as calculating PCA projections or f-statistics) 
due to its robustness with respect to genotyping assumptions. 
To produce such files, you will have to generate pseudo-haploid genotype calls 
from the raw .*bam* file. A common tool used for this purpose is 
*pileupCaller* distributed via the software `sequenceTools <https://github.com/stschiff/sequenceTools>`_ 
Ideally, this pseudo-haploid genotype calling is done in a way adjusted to the data generation process (e.g. clipping a certain amount of terminal base-pairs depending on the form of UDG treatment).

With the *roh.csv* and *roh_full.csv* tables, what is the difference between the *length*/*lengthM*, *start*/*startM* and  *end*/*endM* columns?  
**************************************************************************************************************************************************

The columns without the *M* report the index of SNPs covered in the target. One can use them for quality control (e.g. how many SNPs are covered within an identified ROH) and to quickly identify ROH regions by index in the posterior output (which also only outputs data for covered SNPs). The key result column are the ones with the `M`, which is the length of each ROH, their start and their end position in Morgan.


Why is the value of *lengthM* sometimes multiplied per 100? 
************************************************************

Some of the visualizations and in the full individual table depict the data in centimorgan. Thus, the Morgan value is multiplied with 100 to transform from Morgan to centimorgan.


How can I implement hapROH for another set of SNPs, such as all 1000Genome SNPs?
***********************************************************************************

HapROH requires hdf5 reference files for each chromosome. These hdf5 files have to include the phased genotype data as well as a genetic map for each SNP (as a matching array in the field *variants/MAP*). One can get such files from a .*vcf* file by using the scikit allele function `vcf_to_hdf5`, and adding the genetic map field into the resulting .hdf5. For the latter, one can use the function :meth:`hapsburg.PackagesSupport.h5_python.h5_functions.merge_in_ld_map`.

Example code for preparing the 1240K reference panel including cleaning SNPs and merging in the linkage map can be found `here <https://github.com/hringbauer/hapROH/blob/master/Notebooks/PrepareData/prepare_1000genomes_1240K.ipynb>`_.

For 1000 Genome SNPs, a reference panel (which has not been extensively tested) is available upon request.


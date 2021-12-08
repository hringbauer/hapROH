Introduction
============
The ``hapROH`` package identifies runs of homozygosity (ROH) in ancient and present-day DNA by using a panel of reference haplotypes. This package contains functions and wrappers to call ROH and functions for downstream analysis and visualization.

For downward compatibility, the package uses ``hapsburg`` as module name. After installation you can import Python functions via *from hapsburg.XX import YY*.

Scope of the Method
**********************

Standard parameters are tuned for human 1240K capture data (ca. 1.2 million SNPs used widely in human aDNA analysis) and using 1000 Genome haplotypes as reference. The software is tested on a wide range of test applications, both 1240K data and also whole genome sequencing data downsampled to 1240K SNPs. Successful cases include 45k year old Ust Ishim man, and a wide range of American, Eurasian and Oceanian ancient DNA, showing that the method generally works for split times of reference panel and target up to a few 10k years, which includes all out-of-Africa populations (Attention: Neanderthals and Denisovans do not fall into that range, additionally some Subsaharan hunter gatherer test cases did not give satisfactory results).

Currently, hapROH works on data for 1240K SNPs and in unpacked or packed `eigenstrat` format (which is widely used in human ancient DNA). The software assumes pseudo-haploid or diploid genotype data (the mode can be set, by default it is pseudo-haploid). The recommended coverage range is 400,000 or more 1240K SNPs covered at least once (i.e. at least ca. 0.3x coverage).

If you have whole genome data available, you have to downsample an create eigenstrat files for biallelic 1240k SNPs first.

In case you are planning applications to other kind of SNP or bigger SNP sets, or even other organisms, the method parameters have to be adjusted (the default parameters are specifically optimized for human 1240K data). You can mirror our procedure to find good parameters (described in the publication), and if you contact me for assistance - I am happy to share my own experience.


Citation
**********

If you use the software for a scientific publication and want to cite it, you can use:
https://www.biorxiv.org/content/10.1101/2020.05.31.126912v1

Contact
**********

If you have bug reports, suggestions or general comments, please contact me. I am happy to hear from you. Bug reports and user suggestions will help me to improve this software - so please do not hesitate to reach out!

harald_ringbauer AT eva mpg de
(fill in blanks with dots and AT with @)

Acknowledgments
**********

Big thank you to my co-authors Matthias Steinrücken and John Novembre. The project profited immensely from Matthias' deep knowledge about HMMs and from John's extensive experience in developing population genetics software. Countless discussions with both have been key for moving forward this project. Another big thanks goes to Nick Patterson, who informed me about the benefits of working with rescaled HMMs - substantially improving the runtime of hapROH. 

I want to acknowledge users who find and report software bugs (Mélanie Pruvost, Ke Wang, Ruoyun Hui) and all users who reached out with general questions and requests (Rosa Fregel, Federico Sanchez). This feedback has helped to remove errors in the program and to improve its usability. Many thanks!



Author:
Harald Ringbauer, 2021
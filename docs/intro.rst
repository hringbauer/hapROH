Overview
============
The ``hapROH`` software package contains two primary modules and a range of functions to visualize the results.

#. ``hapROH`` identifies runs of homozygosity (ROH) in ancient and present-day DNA by using a panel of reference haplotypes. As input, it uses genotype files in eigentstrat format. This module contains functions and wrappers to call ROH and functions for downstream analysis and visualization. For downward compatibility, this software uses ``hapsburg`` as module name.

#. ``hapCON`` estimates contamination in aDNA data of male individuals by using a panel of reference haplotypes. It works directly from BAM file or from samtools mpileup output. 


Citations
**********

If you use the software for a scientific publication and want to cite it, you can use:

#. https://doi.org/10.1038/s41467-021-25289-w for ``hapROH``
#. https://doi.org/10.1101/2021.12.20.473429 for ``hapCON``


Contact
**********

If you have bug reports, suggestions or any general comments, please contact us. Please do not hesitate to reach out, we are happy to hear from you! Bug reports and user suggestions will help us to improve this software.

harald_ringbauer AT eva mpg de

yilei_huang AT eva mpg de

(fill in AT with @ and other blanks with dots)

Acknowledgments
*****************

Big thank you to the two original co-authors of hapROH, Matthias Steinrücken and John Novembre. This project and follow-ups profited immensely from Matthias' deep knowledge about HMMs and from John's extensive experience in developing population genetics software. Countless discussions with both have been key for moving forward this project. Another big thanks goes to Nick Patterson, who informed me about the benefits of working with rescaled HMMs - substantially improving the runtime of hapROH. 

We also want to acknowledge everyone who found and reported software bugs (including Mélanie Pruvost, Ke Wang, Ruoyun Hui) and all users who reached out with general questions and requests (including Rosa Fregel, Federico Sanchez). This feedback has helped to remove errors in the program and to improve its usability. Many thanks!


Authors:
Harald Ringbauer, Yilei Huang, 2021

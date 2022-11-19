## 0.62 November 19th, 2022, Harald Ringbauer, Yilei Huang
- Fixed a bug with spawning multiprocessing that prevented hapCon from running.

## 0.61 November 18th, 2022, Harald Ringbauer, Yilei Huang 
[CONTAINS CRITICAL BUG]
- Fixed bug in gap-merging where the bp position was not merged.
- Fixed a hapcon bug caused by a "+" in pileups also denoting deletion.

## 0.60, August 29th 2022, Harald Ringbauer, Yilei Huang
- Included hapCon+ROH functionality from Yilei's branch
- Added downsample parameter to downsample read count data to this specified average coverage

## 0.53, July 1st 2022, Harald Ringbauer
- Merged in Yeilei's branch, that includes Unix-style command line functionality for hapCon

## 0.52, June 29th 2022, Harald Ringbauer
- Updated the license (so that license file is automatically available) and package meta information, added setup.cfg for that.

## 0.51a0, June 8th 2022, Harald Ringbauer
- Minor bug fix in h5_functions (Import now works again - fixing the issue with
  with the incorrectly formatted `major_foc.0f`

## 0.5a, December 1st 2022, Yilei Huang
This is the addition of hapCON - estimating contamination in aDNA with haplotype copying!
- added hapCON: a new method for estimating male X chromosome contamination
- added new parameters for hapsb_ind and hapsb_chrom to perform automatic merging of short ROH blocks into a long block

## 0.3a4, June 14th 2021, Harald Ringbauer
- Fixed minor bug leading to posterior for last locus on target chromosome not being normalized probably.

## 0.3a3, April 26th 2021, Harald Ringbauer
- Updates to the README, including a new FAQ section.

## 0.3a2, March 19th 2021, Harald Ringbauer
Minor improvements and quality of life updates:
- The hapsb_chrom now plots all of its function variables. This will be printed into the logfile - the idea is to have better reproducibiliyt.
- Implemented some runtime updates for the read count model. There is now an Einstein sum done when summing over the latent genotypes. This prevents extra space being used. Moreover, readcount data are now int8 by default. Attention: This limits the maximal read count for each locus to 127!
- Fixed a bug in the diploid_gt emission mode where an undefined e_mat0 was used.

## 0.3a1, December 14th 2020, Harald Ringbauer
Major runtime update:
- The HMM is now rescaling every SNP instead of working in SNP space. This makes the HMM 8x, and the overall runtime ~60% faster. The main bottleneck now is loading the reference data. 
- For backward compatibility, the posterior plotting function automatically detects log/non-log space.
- Implemented that all reference SNPs are loaded into RAM (fast for HDF5), but that only the SNPs that are covered are extracted, in a one-step downsampling (instead of two step downsampling before.)

## 0.2a3, November 18th 2020, Harald Ringbauer
- Introduced a function to plot a histogram of multiple indiviudals and expectation of multiple indivdiuals (plot_pde_individual can now take a list of individuals as input). 
- Fixed typo in this function name: It is now plot_pde_individual (previously: plot_pde_indivdiual).

## 0.2a2, November 18th 2020, Harald Ringbauer
- Hot fix to deal with issue with new h5py-3.1.0. Now hdf5 char strings from the REF/ALT column are loaded as bytestrings that don't match the strings from eigenstrat. Implemented a conversion to string.

## 0.2a1, October 19th 2020, Harald Ringbauer
- New additional output (after several requests from users): Give out Physical Positions in coordinates of GRCh37 in hapROH output. 
  These are additional columns, so they do not change any expected behavior. Upgraded version to have clean cut on output change.
- Internal: Added tests for some plotting functions to check expected behavior.

## 0.1a9, October 1st 2020, Harald Ringbauer
- Add Requirement for psutil package in installation (needed to print memory usage)
- Fixed bug with a lost string "random" in loadeigenstrat.

## 0.1a8, September 1st 2020, Harald Ringbauer
- Added function in pp_individual_roh_csvs to post-process single hapROH ROH data-frame. It creates combined dataframes with every single ROH row.
- Updated README.md, reflecting that additional vignettes (for calling IBD on the X as well as pop size estimation) have been added.

## 0.1a7, August 24th 2020, Harald Ringbauer
- Added likelihood profile inference for Ne estimator, including for 95% CI intervalls based on 1.92 LL units down the maximum.
- Added support for different chromosome lengths for Ne estimator

## 0.1a6, July 23th 2020, Harald Ringbauer
- Fixed bug that for big eigenstrat sizes the automatic eigenstrat is_binary checker got stuck.
- Minor quality of life improvements for X data. Improved and more consistent output (using chr23 instead of chrX)

## 0.1a5, July 22th 2020, Harald Ringbauer
- Introduced mode `Eigenstrat` for preprocessing. It automatically detects whether an Eigenstrat is packed or unpacked.
  The package loadEigenstrat has a packed=-1 in the function to allow for automated choice.
- Introduced preprocessing mode `EigenstratX` and postprocessing mode `IBD_X` that allows to call IBD on the X between two males, using pseudo-haploid eigenstrat and a X Chromosome reference panel as input. Also introduced a function in hapsburg_run (`hapsb_chromXs`) that can run multiple X chromsosomes, with multi-processing. Also a special post-processing function designed for X Chromosome data. (pp_X_roh in pp_individual_roh_csvs)

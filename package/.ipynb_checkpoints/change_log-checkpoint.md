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

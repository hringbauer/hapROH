## 0.1a5, July 10th 2020, Harald Ringbauer

- Introduced mode `Eigenstrat` for preprocessing. It automatically detects whether an Eigenstrat is packed or unpacked.
  The package loadEigenstrat has a packed=-1 in the function to allow for automated choice.
- Introduced preprocessing mode `EigenstratX` and postprocessing mode `IBD_X` that allows to call IBD on the X between two males, using pseudo-haploid eigenstrat and a X Chromosome reference panel as input. Also introduced a function in hapsburg_run (`hapsb_chromXs`) that can run multiple X chromsosomes, with multi-processing. Also a special post-processing function designed for X Chromosome data. (pp_X_roh in pp_individual_roh_csvs)
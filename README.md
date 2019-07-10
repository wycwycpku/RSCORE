# RSCORE
RSCORE is an R package to explore single-cell RNA-seq data with the view of molecular networks.

To run RSCORE, you need to install some extra dependencies:
```
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
## Required
BiocManager::install(c("Seurat", "AUCell", "STRINGdb", "propr", "coop", "network"))
## To support paralell execution (not available in Windows)
BiocManager::install(c("doMC", "doRNG"))
## To get marker genes quickly
BiocManager::install("mahmoudibrahim/genesorteR") 
```

Then you can install RSCORE:
```
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("wycwycpku/RSCORE")
```
if you want to build vignettes, you need to add an option
```
devtools::install_github("wycwycpku/RSCORE", build_opts = NULL)
```
More details and tutorials can be found at:
[Introduction and Tutorials](https://github.com/wycwycpku/RSCORE/blob/master/vignettes/RSCORE_Tutorials.Rmd)
	
# RSCORE
RSCORE is an R package to explore single-cell RNA-seq data with the view of molecular networks. The method is outlined in our manuscript
**Enhancing single-cell cellular state inference by incorporating molecular network features**
The preprint version is available on [bioRxiv](https://doi.org/10.1101/699959). 

To run RSCORE, you need to install some extra dependencies:
```
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
## Required
BiocManager::install(c("Seurat", "AUCell", "STRINGdb", "propr", "coop", "network"))
## To support parallel execution
BiocManager::install(c("doMC", "doRNG"))
## To get marker genes quickly
BiocManager::install("mahmoudibrahim/genesorteR") 
## To do GO enrichment
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
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
[Introduction and Tutorials](https://github.com/wycwycpku/RSCORE/blob/master/vignettes/RSCORE_Tutorials.pdf)
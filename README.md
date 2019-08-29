# RSCORE
RSCORE is an R package to explore single-cell RNA-seq data with the view of molecular networks. The method is outlined in our manuscript 

**Enhancing single-cell cellular state inference by incorporating molecular network features**

*Ji Dong, Peijie Zhou, Yichong Wu, Wendong Wang, Yidong Chen, Xin Zhou, Haoling Xie, Jiansen Lu, Xiannian Zhang, Lu Wen, Wei Fu, Tiejun Li, Fuchou Tang*.

The preprint version is available on [bioRxiv](https://doi.org/10.1101/699959). 

![The workflow of RSCORE](https://github.com/wycwycpku/RSCORE/blob/master/images/Figure_1.jpg) 
## Install
To run RSCORE, you need to install some extra dependencies:
```
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
## Required
BiocManager::install(c("Seurat", "AUCell", "STRINGdb", "propr", "coop", "network","intergraph"))
## To support parallel execution
BiocManager::install("doRNG")
install.packages("doMC", repos="http://R-Forge.R-project.org")
install.packages("doParallel")
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

## Tutorial
For more details and basic usage see following tutorials
1.	[Single cells from different cell lines](https://github.com/wycwycpku/RSCORE/blob/master/vignettes/RSCORE_Tutorials.pdf)
2.	[Human embryonic cell altas (HECA)](https://github.com/zorrodong/HECA)

## News
2019/07/23:
-	Add GO enrichment function.
-	Correct some bugs.

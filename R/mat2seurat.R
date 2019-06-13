#' processing function
#'
#' transfer raw matrix data to Seurat object
#'
#' @param Data gene*cell matrix
#' @param min_nFeature Include cells where at least this many features are detected
#' @param min_nCount Include cells where at least this many counts are detected
#' @param features Vector of features to use as variable features.
#' @param selection.method How to choose top variable features.
#' @param nfeatures Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'
#' @param mean.cutoff A two-length numeric vector with low- and high-cutoffs for feature means
#' @param dispersion.cutoff A two-length numeric vector with low- and high-cutoffs for feature dispersions
#' @param names.filed
#' @param names.delim

#' @return Seurat object
#' @export
#'
#' @examples
mat2seurat <- function( Data, min_nFeature = 1000, min_nCount = 10000,
                        features = NULL, selection.method = c('mean.var.plot', 'vst', 'dispersion'),
                        nfeatures = NULL, mean.cutoff = c(0.1,Inf), dispersion.cutoff = c(0.1,Inf),
                        names.filed = 1, names.delim = '_' )
{
  Data <- CreateSeuratObject(counts = as.matrix(Data), project = 'Data', assay = 'RNA',
                             min.cells = 3, min.features = min_nFeature)
  eval(parse(text= paste('Data <- subset(x = Data, subset =  nCount_RNA >', min_nCount,')')))
  Data <- NormalizeData(object = Data, scale.factor = 100000)

  ## feature selection methods can be adjusted
  selection.method <- selection.method[1]
  if(!is.null(features)){
    VariableFeatures(Data) <- features
  }else if(!is.null(nfeatures)){
    Data <- FindVariableFeatures(Data, selection.method = selection.method, nfeatures = nfeatures)
  }else
    Data <- FindVariableFeatures(Data, selection.method = selection.method, mean.cutoff = mean.cutoff,
                                 dispersion.cutoff = dispersion.cutoff)
  Data <- ScaleData(object = Data)

  return(Data)
}

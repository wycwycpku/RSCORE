#' Find expression markers for all identity classes
#'
#' Find expression markers for all identity classes
#'
#' @param object Seurat object
#' @param assay assay of the object, default is the default assay
#' @param binarizeMethod Either "median" (default) or "naive" or a numeric cutoff
#' @param FoldChange default is 2
#' @param p.adj default is 0.05
#' @param cores default is 1
#'
#' @return
#' @export
#' @import genesorteR
#'
#' @examples
Find_Markers <- function(object, assay = NULL, binarizeMethod = "median", FoldChange = 2,
                         p.adj = 0.05, cores = 1){
  if(is.null(assay)){
    assay <- DefaultAssay(object)
  }
  expr_mtx <- GetAssayData(object = object, assay = assay)
  gs <- sortGenes(expr_mtx, Idents(object), binarizeMethod = binarizeMethod, cores = cores)
  pp <- getPValues(gs, cores = cores)

  gene_score <- gs$specScore
  gene_pvalue <- pp$pval
  gene_adj.p <- pp$adjpval

  Markers <- data.frame()
  for(i in 1:ncol(gene_adj.p)){
    temp_gene <- rownames(gene_adj.p)[gene_adj.p[,i] < p.adj]

    if(length(temp_gene) != 0){
      temp_cluster <- colnames(gene_adj.p)[i]

      temp_cluster_cells <- WhichCells(object = object, idents = temp_cluster)
      temp_cluster_mtx <- expr_mtx[temp_gene,temp_cluster_cells]

      temp_other_cells <- WhichCells(object = object, idents = colnames(gene_adj.p)[colnames(gene_adj.p)!=temp_cluster])
      temp_other_mtx <- expr_mtx[temp_gene,temp_other_cells]

      if(length(temp_gene)>1){
        temp_fc <- Matrix::rowMeans(temp_cluster_mtx)/Matrix::rowMeans(temp_other_mtx)
        temp_gene_filter <- temp_gene[temp_fc >= FoldChange]
        temp_fc_filter <- temp_fc[temp_fc >= FoldChange]
        temp_gene_sort <- names(sort(gene_score[temp_gene_filter,i],decreasing = T))
      }else{
        temp_fc <- mean(temp_cluster_mtx)/mean(temp_other_mtx)
        temp_gene_sort <- temp_gene
        temp_fc_filter <- temp_fc[temp_fc >= FoldChange]
      }

      if(!is.null(temp_gene_sort)){
        temp_matrix <- data.frame(Marker = temp_gene_sort,
                                  Cluster = colnames(gene_adj.p)[i],
                                  FoldChange = temp_fc_filter,
                                  P.value = gene_pvalue[temp_gene_sort,i],
                                  P.adj = gene_adj.p[temp_gene_sort,i],
                                  Gene.Score = gene_score[temp_gene_sort,i])
        Markers <- rbind(Markers,temp_matrix)
      }
    }
  }
  rownames(Markers) <- 1:nrow(Markers)
  Markers$Marker <- as.character(Markers$Marker)
  Markers$Cluster <- as.character(Markers$Cluster)

  Markers_list <- list(Markers=Markers,GeneSort=gs)

  return(Markers_list)
}


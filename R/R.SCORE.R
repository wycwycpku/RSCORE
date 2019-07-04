#' R.SCORE
#'
#' explore single-cell RNA-seq data with the view of molecular networks.
#'
#' @param Data gene*cell matrix or a seurat object.
#' @param PPI protein-protein network, you can provide a matrix if you have PPI data,
#'            or you can let PPI = 'Biogrid' or 'String', the system will download corresponding PPI data.
#' @param species NCBI taxonomy identifier of the species to query for homologs, default is 9606(Homosapiens),
#'                only needed when PPI = 'Biogrid' or 'String'.
#' @param score_threshold threshold on the score of the interaction, default is 600, only needed when PPI = 'String'
#' @param metric A character string. The metric used to calculate association.
#'               Choose one of 'cor', 'rho', 'phi', 'phs'
#' @param module_min Integer number. The minimum size of a module to be accepted.
#' @param max_step Integer number. The maximum step run in the community detect.
#' @param AUCRank Integer number. The number of highly-expressed genes to include when computing AUCell.
#' @param nCores Number of cores to use for computation.
#'
#' @return a seurat object
#' @export
#' @import Seurat AUCell igraph propr
#' @importFrom coop pcor
#'
#' @examples
R.SCORE <- function(Data, PPI = 'Biogrid', species = 9606, score_threshold = 600,
                  metric = c('cor','rho','phi','phs'),
                  module_min = 3, max_step = 10, AUCRank = 400, nCores = 1)
{
  ## data processing
  if(class(Data) == 'matrix'){
    Data <- mat2seurat(Data)
  }else if(class(Data) != 'Seurat'){
    stop("Wrong data class")
  }

  ## get PPI network
  if(class(PPI) == 'matrix'){
    hs_network_matrix <- PPI
  }else if(PPI == 'String'){
    hs_network_matrix <- getPPI_String(Data, species = species, score_threshold = score_threshold)
  }else if(PPI == 'Biogrid'){
    hs_network_matrix <- getPPI_Biogrid(Data, species = species)
  }else{
    stop("Provided 'PPI' not recognized")
  }
  rm(PPI)

  Data_expr <- as.matrix(GetAssayData(Data))
  Data_expr_scale <- as.matrix(GetAssayData(Data, slot = 'scale.data'))
  var_genes <- VariableFeatures(Data)
  Data_genes <- var_genes[var_genes %in% as.character(rownames(hs_network_matrix))]

  ## update the expression matrix
  Data_expr <- Data_expr[Data_genes,]
  Data_expr_scale <- Data_expr_scale[Data_genes,]

  ## calculate the correlation coefficient or other metric
  metric <- metric[1]
  if(metric == 'cor')
  {
    Data_metric <- coop::tpcor(Data_expr)
    Data_metric[Data_metric < 0] <- 0
  }else{
    Data_expr_counts <- as.matrix(GetAssayData(Data, slot = 'counts'))
    Data_expr_counts <- Data_expr_counts[Data_genes,]
    if(metric == 'rho')
    {
      Data_metric <- getMatrix(propr(t(Data_expr_counts),metric = 'rho'))
    }else if(metric == 'phi')
    {
      Data_metric <- getMatrix(propr(t(Data_expr_counts),metric = 'phi'))
    }else if(metric == 'phs')
    {
      Data_metric <- getMatrix(propr(t(Data_expr_counts),metric = 'phs'))
    }else{
      stop("Provided 'metric' not recognized.")
    }
    rm(Data_expr_counts)
  }

  ## construct weighted network
  hs_network_final <- hs_network_matrix[Data_genes,Data_genes]
  Data_network_final <- Data_metric * hs_network_final
  filtered_row <- rowSums(as.matrix(Data_network_final)) > 0
  Data_network_final <- Data_network_final[filtered_row,filtered_row]
  rm(Data_metric)
  rm(hs_network_matrix)
  rm(hs_network_final)

  network_trim <- igraph::graph_from_adjacency_matrix(Data_network_final,weighted = TRUE)

  ## network decomposition
  gene_sets_all <- list()
  run_label = c()
  for (rw_step in 1:max_step){
    set.seed(123)
    network_cluster <- igraph::walktrap.community(network_trim,steps =rw_step)
    gene_sets_all <- c(gene_sets_all,communities(network_cluster))
    run_label <- c(run_label,rep(rw_step,length(communities(network_cluster))))
  }

  temp_len <- c()
  geneSets <- list()
  k = 0;
  for (i in 1:length(gene_sets_all)){
    temp <- as.character(unlist(gene_sets_all[i]))
    if(length(temp)>=module_min){
      k <- k+1;
      temp_len[k] <- length(temp)
      geneSets[[paste('Set',run_label[i],k,temp_len[k],sep = ".")]] <- temp
    }
  }
  cat('module num: ',length(geneSets),'\n')

  ## calculate AUC scores
  cells_rankings <- AUCell_buildRankings(Data_expr_scale, nCores = nCores)
  rm(Data_expr)
  rm(Data_expr_scale)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=AUCRank, nCores = nCores)
  cells_AUC_matrix <- getAUC(cells_AUC)

  ## creater net-based seurat object
  Data[["Net"]] <- CreateAssayObject(counts = as.matrix(cells_AUC_matrix))
  Misc(Data, slot = 'geneSets') <- geneSets
  Misc(Data, slot = 'Data_network') <- Data_network_final
  DefaultAssay(Data) <- "Net"
  Data <- ScaleData(Data, assay = "Net")

  return(Data)
}

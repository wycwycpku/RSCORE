#' Title
#'
#' @param object an Seurat object which contains the expression matrix
#' @param version String version
#' @param species NCBI taxonomy identifier of the species to query for homologs, default is 9606(Homosapiens)
#' @param score_threshold threshold on the score of the interaction, default is 600
#'
#' @return PPI matrix
#' @import STRINGdb
#' @export
#'
#' @examples
getPPI_String <- function(object, version = '10', species = 9606, score_threshold = 600)
{
  string_db <- STRINGdb$new(version = version, species = species, score_threshold = score_threshold, input_directory = '')
  Data_gene <- as.data.frame(rownames(object))
  colnames(Data_gene)[1] <- 'geneid'
  Data_gene$geneid_backup <- Data_gene$geneid
  Data_gene <- string_db$map(Data_gene, 'geneid', removeUnmappedRows = T)
  Data_gene <- as.data.frame(unique(Data_gene, by='STRING_id'))
  stringid <- Data_gene[,'STRING_id']

  hs_network <- string_db$get_graph()
  hs_network_matrix <- as_adj(hs_network)
  rowid <- stringid[stringid %in% rownames(hs_network_matrix)]
  genename <- Data_gene[stringid %in% rownames(hs_network_matrix),'geneid_backup']
  hs_network_matrix <- hs_network_matrix[rowid,rowid]
  rownames(hs_network_matrix) <- genename
  colnames(hs_network_matrix) <- genename
  saveRDS(hs_network_matrix, 'hs_network_matrix_String10.Rda')

  return(hs_network_matrix)
}

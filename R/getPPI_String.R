#' get PPI from String
#'
#' get PPI from String
#'
#' @param object an Seurat object which contains the expression matrix
#' @param version String version
#' @param species NCBI taxonomy identifier of the species to query for homologs, default is 9606(Homosapiens)
#' @param score_threshold threshold on the score of the interaction, default is 600
#' @param save whether save the result as a Rdata file, default is FALSE
#'
#' @return PPI matrix
#' @import STRINGdb data.table
#' @export
#'
#' @examples
getPPI_String <- function(object, version = '10', species = 9606, score_threshold = 600, save = FALSE)
{
  string_db <- STRINGdb$new(version = version, species = species, score_threshold = score_threshold, input_directory = '')

  Data_gene <- as.data.frame(rownames(object))
  colnames(Data_gene)[1] <- 'geneid'
  Data_gene$geneid_backup <- Data_gene$geneid
  Data_gene <- string_db$map(Data_gene, 'geneid', removeUnmappedRows = T)
  Data_gene <- as.data.frame(unique(as.data.table(Data_gene), by='STRING_id'))
  Data_gene <- as.data.frame(unique(as.data.table(Data_gene), by='geneid_backup'))
  stringid <- Data_gene[,'STRING_id']

  net <- string_db$get_graph()
  net_adj <- as_adj(net)

  rowid <- stringid[stringid %in% rownames(net_adj)]
  genename <- Data_gene[stringid %in% rownames(net_adj),'geneid_backup']
  net_adj <- net_adj[rowid,rowid]
  rownames(net_adj) <- genename
  colnames(net_adj) <- genename

  if(save){
    saveRDS(net_adj, paste(species,'_ppi_matrix_STRING-',version,'.Rda',sep=""))
  }

  return(net_adj)
}

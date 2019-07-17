#' get PPI from Biogrid
#'
#' get PPI from Biogrid
#'
#' @param object
#' @param species NCBI taxonomy identifier of the species to query for homologs, default is 9606(Homosapiens)
#' @param biogrid_version Version of BioGRID, such as '3.5.174'
#'
#' @return PPI matrix
#' @export
#'
#' @examples
getPPI_Biogrid <- function(object = NULL, biogrid_version='3.5.174', species = 9606)
{
  dbFiles <- paste('https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-',biogrid_version,'/BIOGRID-ALL-',biogrid_version,'.tab2.zip',sep="")
  if(!file.exists(sub(pattern = 'zip', replacement = 'txt', x = basename(dbFiles))))
  {
    if(!file.exists(basename(dbFiles)))
      download.file(dbFiles, destfile = basename(dbFiles))
    unzip(basename(dbFiles),exdir = '.')
  }
  PPI <- read.csv(sub(pattern = 'zip', replacement = 'txt', x = basename(dbFiles)), header = T, sep = '\t')
  species_coincide <- PPI[,'Organism.Interactor.A'] == species & PPI[,'Organism.Interactor.B'] == species & PPI[,'Experimental.System.Type'] == 'physical'
  PPI = PPI[species_coincide,]

  if(!is.null(object))
  {
    gene_data <- rownames(object)
    PPI <- PPI[which(is.element(PPI[,'Official.Symbol.Interactor.A'], gene_data)),]
    PPI <- PPI[which(is.element(PPI[,'Official.Symbol.Interactor.B'], gene_data)),]
  }
  gene_name <- union(PPI[,'Official.Symbol.Interactor.A'],PPI[,'Official.Symbol.Interactor.B'])
  nodes <- gene_name
  links <- PPI[,c('Official.Symbol.Interactor.A','Official.Symbol.Interactor.B')]
  net <- graph_from_data_frame(d = links,vertices = nodes,directed = FALSE)
  net <- igraph::simplify(net)
  saveRDS(as_adj(net),paste('hs_ppi_matrix_BioGRID-',biogrid_version,'.Rda',sep=""))
  file.remove(paste('BIOGRID-ALL-',biogrid_version,'.tab2.txt',sep=""))
  file.remove(paste('BIOGRID-ALL-',biogrid_version,'.tab2.zip',sep=""))
  return(as_adj(net))
}

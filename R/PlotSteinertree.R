#' Title
#'
#' @param object
#' @param ident which class you want to show, default: all
#' @param weighted if it is TRUE, the network will be weighted, default TRUE
#' @param onlysteiner if it is TRUE, will only show the steiner tree, default TRUE
#'
#' @return
#' @export
#'
#' @examples
PlotSteinertree <- function(object, ident = NULL, weighted = T, onlysteiner = T)
{
  Data_corr <- object@misc$Data_net
  Data_dist <- 1/(Data_corr + 1e-3)

  #weighted net
  network_trim <- graph_from_adjacency_matrix(Data_dist, mode = 'undirected', weighted = T)
  connected <- components(network_trim)
  genes_in_connected_set <- names(connected$membership[connected$membership == 1])

  if(is.null(ident)){
    RNA.markers <- FindAllMarkers(object, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
    Net.markers <- FindAllMarkers(object, assay = "Net", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
    type <- levels(Idents(object))
  }else{
    RNA.markers <- FindMarkers(object, assay = "RNA", ident.1 = ident, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
    Net.markers <- FindMarkers(object, assay = "Net", ident.1 = ident, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
    type <- ident
  }
  for(i in type){
    #Select genes
    if(is.null(ident)){
      genes_in_select_markers <- rownames(RNA.markers[RNA.markers$cluster %in% i  & RNA.markers$p_val_adj < 0.01,])
      modules_select <- rownames(Net.markers[Net.markers$cluster %in% i,])
    }else{
      genes_in_select_markers <- rownames(RNA.markers[RNA.markers$p_val_adj < 0.01,])
      modules_select <- rownames(Net.markers)
    }
    marker_genes <- genes_in_select_markers[genes_in_select_markers %in% genes_in_connected_set]
    modules_name_all <- object@misc$geneSets
    genes_in_select_modules <- unique(rownames(table(unlist(modules_name_all[modules_select]))))
    module_genes <- genes_in_select_modules[genes_in_select_modules %in% genes_in_connected_set]

    both_genes <- marker_genes[marker_genes %in% module_genes]
    terminals <- unique(c(marker_genes, module_genes))

    steiner_tree_sp <- steinertree(terminals = terminals, graph = network_trim, weighted = weighted)

    whichfig <- onlysteiner + 1
    ###Color
    V(steiner_tree_sp[[whichfig]])$color <- 'white'
    V(steiner_tree_sp[[whichfig]])[marker_genes]$color <- 'red'
    V(steiner_tree_sp[[whichfig]])[module_genes]$color <- 'green'
    V(steiner_tree_sp[[whichfig]])[both_genes]$color <- 'yellow'

    set.seed(0)
    plot(steiner_tree_sp[[whichfig]], vertex.size=4, main=i, edge.curved=0.3,
         vertex.label.cex=.5,vertex.label.dist=1)
  }
}

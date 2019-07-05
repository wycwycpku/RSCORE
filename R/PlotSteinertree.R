#' Plot steiner tree
#'
#' show the steiner tree of given clusters
#'
#' @param object
#' @param ident which class you want to show, default: all
#' @param weighted if it is TRUE, the network will be weighted, default TRUE
#' @param onlysteiner if it is TRUE, will only show the steiner tree, default TRUE
#' @param method how to find steinertree, default is 'sp' (Shortest Path Based Approximation)
#'               'kb' (Kruskal's minimum spanning tree algorithm) also can be choosen. Be careful, 'kb' is very slow.
#' @param nodes the nodes you want to labeled in the Steiner tree
#' @param typeof_node_size 'pagerank' or 'degree', former is default
#' @param only_label_terminal default is F
#' @param label_num number of nodes you want to label, default is 20
#' @return
#' @export
#'
#' @examples
PlotSteinertree <- function(object, ident = NULL, nodes = NULL, method = 'sp', weighted = T, onlysteiner = T,
                            typeof_node_size = c('pagerank','degree'), only_label_terminal = F, label_num = 20)
{
  Data_corr <- object@misc$Data_net
  network_trim <- igraph::graph_from_adjacency_matrix(Data_corr, mode = 'undirected', weighted = T)
  Data_dist <- 1/(Data_corr + 1e-3)

  #weighted net
  network_trim <- igraph::graph_from_adjacency_matrix(Data_dist, mode = 'undirected', weighted = T)
  connected <- igraph::components(network_trim)
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
  plots <- list()
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

    both_genes <- intersect(marker_genes, module_genes)
    terminals <- union(marker_genes, module_genes)
    #terminals <- union(terminals, nodes)

    cat('calculate tree\n')
    steiner_tree_sp <- steinertree(terminals = terminals, graph = network_trim, method= method, weighted = weighted)

    whichfig <- onlysteiner + 1
    plots[[as.character(i)]] <- get_steiner_plot(steiner_tree_sp[[whichfig]], label_num = label_num, name = i, typeof_node_size = typeof_node_size,
                                                 terminals = terminals, marker_genes = marker_genes, module_genes = module_genes,
                                                 nodes = nodes, only_label_terminal = only_label_terminal)

  }
  return(plots)
}

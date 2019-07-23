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
#' @param nodes the nodes you want to labeled in the Steiner tree, if it is NULL (default), see 'label_num'
#' @param typeof_node_size 'pagerank' or 'degree', former is default
#' @param only_label_terminal default is F
#' @param geneset1 If ident == NULL, then you should provide two set of genes to compare
#' @param geneset2 same to geneset1
#' @param label_num If nodes == NULL, the number of nodes you want to label, default is 20
#'
#' @return
#' @export
#'
#' @examples
PlotSteinertree <- function(object, ident = NULL, geneset1 = NULL, geneset2 = NULL, nodes = NULL, method = 'sp',
                            weighted = T, onlysteiner = T, typeof_node_size = c('pagerank','degree'),
                            only_label_terminal = F, label_num = 20)
{
  Data_corr <- as.matrix(object@misc$Data_net)
  network_trim <- igraph::graph_from_adjacency_matrix(Data_corr, mode = 'undirected', weighted = T)
  Data_dist <- 1/(Data_corr + 1e-3)

  #weighted net
  network_trim <- igraph::graph_from_adjacency_matrix(Data_dist, mode = 'undirected', weighted = T)
  connected <- igraph::components(network_trim)
  genes_in_connected_set <- names(connected$membership[connected$membership == 1])

  if(is.null(ident)){
    if(is.null(geneset1) || is.null(geneset2)){
      stop('You should provide the ident of one cluster or two set of genes.\n')
    }else{
      gs1_genes <- intersect(geneset1, genes_in_connected_set)
      gs2_genes <- intersect(geneset2, genes_in_connected_set)
      if(length(gs1_genes) == 0){
        stop('There is no gene of geneset1 in the network!\n')
      }
      if(length(gs2_genes) == 0){
        stop('There is no gene of geneset2 in the network!\n')
      }
    }
    name <- 'Steiner tree'
  }else{
    if(!is.null(geneset1) || !is.null(geneset2)){
      warning('You have provided the ident of cluster, then the geneset is unnecessary.\n')
    }
    SCORE_DEGs_list <- Find_Markers(object = object, assay = 'RNA', FoldChange = 1.5)
    SCORE_DAMs_list <- Find_Markers(object = object, assay = 'Net', FoldChange = 1.5)

    DEGs <- SCORE_DEGs_list$Markers[SCORE_DEGs_list$Markers$Cluster==ident,]$Marker
    DAMs <- SCORE_DAMs_list$Markers[SCORE_DAMs_list$Markers$Cluster==ident,]$Marker
    DAMGs <- unique(rownames(table(unlist(object@misc$geneSets[DAMs]))))
    gs1_genes <- DEGs[DEGs %in% genes_in_connected_set]
    gs2_genes <- DAMGs[DAMGs %in% genes_in_connected_set]
    name <- ident
  }
  both_genes <- intersect(gs1_genes, gs2_genes)
  terminals <- union(gs1_genes, gs2_genes)

  cat('calculate tree\n')
  set.seed(100)
  steiner_tree_sp <- steinertree(terminals = terminals, graph = network_trim, method= method, weighted = weighted)

  whichfig <- onlysteiner + 1
  pl <- get_steiner_plot(steiner_tree_sp[[whichfig]], label_num = label_num, name = name,
                         typeof_node_size = typeof_node_size, terminals = terminals,
                         marker_genes = gs1_genes, module_genes = gs2_genes,
                         nodes = nodes, only_label_terminal = only_label_terminal)
  st_res <- list()
  st_res$tree <- steiner_tree_sp[[whichfig]]
  st_res$plot <- pl
  return(st_res)
}

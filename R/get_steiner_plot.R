#' Title
#'
#' @param ig_obj the Steiner tree, an igraph object
#' @param label_num number of nodes you want to label, default is 20
#' @param name name of the plot
#' @param typeof_node_size 'pagerank' or 'degree', former is default
#' @param only_label_terminal default is F
#' @param terminals terminals of the Steiner tree
#' @param marker_genes marker_genes of the Steiner tree
#' @param module_genes module_genes of the Steiner tree
#' @param nodes the nodes you want to labeled in the Steiner tree
#'
#' @return
#'
#' @import ggrepel ggplot2
#' @importFrom network as.matrix.network.adjacency as.matrix.network.edgelist get.vertex.attribute
#' @examples
get_steiner_plot <- function(ig_obj, label_num = 20, name, typeof_node_size = c('pagerank','degree'), only_label_terminal = T,
                             terminals, marker_genes = NULL, module_genes = NULL, nodes = NULL)
{
  pr_value <- igraph::page_rank(ig_obj, directed = F)
  ig_obj <- igraph::set_vertex_attr(ig_obj, "pagerank", value = scale(pr_value$vector, center = min(pr_value$vector)))
  degrees <- igraph::degree(ig_obj, normalized=FALSE)
  ig_obj <- igraph::set_vertex_attr(ig_obj, "degree", value = degrees)

  max_n <- min(label_num, length(degrees))
  net_obj <- intergraph::asNetwork(ig_obj)
  m <- network::as.matrix.network.adjacency(net_obj) # get sociomatrix
  # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
  # plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))
  # or get it them from Kamada-Kawai's algorithm:
  # plotcord <- data.frame(sna::gplot.layout.kamadakawai(m, NULL))
  plotcord <- data.frame(igraph::layout_nicely(ig_obj))


  colnames(plotcord) <- c("X1","X2")
  edglist <- network::as.matrix.network.edgelist(net_obj)
  edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
  plotcord$vertex.names <- as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
  plotcord$Degree <- network::get.vertex.attribute(net_obj, "degree")
  plotcord$Pagerank <- network::get.vertex.attribute(net_obj, "pagerank")
  plotcord[, "shouldLabel"] <- FALSE
  plotcord[, "Hub"] <- "Steiner genes"
  terminals_bool <- plotcord[,"vertex.names"] %in% terminals
  temp <- plotcord[which(terminals_bool),]
  typeof_node_size <- typeof_node_size[1]
  if (typeof_node_size == 'pagerank') {
    if (only_label_terminal){
      terminals <- temp[order(temp$Pagerank,decreasing = T)[1:max_n],"vertex.names"]
    }else{
      terminals <- plotcord$vertex.names[order(plotcord$Pagerank,decreasing = T)[1:max_n]]
    }
  }else if (typeof_node_size == 'degree') {
    if (only_label_terminal){
      terminals <- temp[order(temp$Degree,decreasing = T)[1:max_n],"vertex.names"]
    }else{
      terminals <- plotcord$vertex.names[order(plotcord$Degree,decreasing = T)[1:max_n]]
    }
  }else{
    stop('Unidentified type of node size!')
  }

  if(!is.null(marker_genes)){
    Both <- intersect(marker_genes, module_genes)
    DEGs <- setdiff(marker_genes, module_genes)
    DAMGs <- setdiff(module_genes, marker_genes)
    Both_bool <- plotcord[,"vertex.names"] %in% Both
    DEGs_bool <- plotcord[,"vertex.names"] %in% DEGs
    DAMGs_bool <- plotcord[,"vertex.names"] %in% DAMGs
    plotcord[which(Both_bool), "Hub"] <- "Both"
    plotcord[which(DEGs_bool), "Hub"] <- "DEGs"
    plotcord[which(DAMGs_bool), "Hub"] <- "DAMGs"
  }

  if(is.null(nodes)){
    sel_vertex <- terminals
  }else{
    cat(sum(nodes %in% plotcord$vertex.names),'/',length(nodes),'of given nodes are in the Steiner tree.\n')
    sel_vertex <- nodes
  }

  colnames(edges) <-  c("X1","Y1","X2","Y2")
  plotcord[which(plotcord[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <- TRUE

  pl <- ggplot(plotcord)  +
    geom_curve(data=edges, aes_(x=~X1, y=~Y1, xend=~X2, yend=~Y2), curvature = 0.3,
               size = 0.5, alpha=0.5, colour="#DDDDDD") +
    geom_label_repel(aes_(x=~X1, y=~X2, label=~vertex.names, color=~Hub),show.legend=F,
                     box.padding=unit(1, "lines"),
                     data=function(x){x[x$shouldLabel, ]}) +
    scale_colour_manual(values=c("Both" = "#E41A1C",
                                 "DAMGs" = "#4DAF4A",
                                 "DEGs" = "#FF7F00",
                                 "Steiner genes" = "#999999")) +
    labs(title=name) +
    ggplot2::theme_bw(base_size = 12, base_family = "") +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = "white",
                                                            colour = NA),
                   panel.border = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  if(typeof_node_size == 'degree'){
    pl <- pl + geom_point(aes_(x=~X1, y=~X2, size=~Degree, alpha=~Degree, color=~Hub))
  }else{
    pl <- pl + geom_point(aes_(x=~X1, y=~X2, size=~Pagerank, alpha=~Pagerank, color=~Hub)) +
      scale_alpha_continuous(breaks=c(1,2),labels=c('low','high')) +
      scale_size_continuous(breaks = c(1,2),labels=c('low','high'))
  }
  return(pl)
}

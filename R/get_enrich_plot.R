#' GO enrichment of Steiner tree
#'
#' GO enrichment of specified group in Steiner tree
#'
#' @param st_res result of cut_steiner_tree (or you can add group information manually)
#' @param group the group of genes you want to do GO enrichment
#'
#' @return a list with enrich_res and plot
#' @import clusterProfiler org.Hs.eg.db cowplot
#' @export
#'
#' @examples
get_enrich_plot <- function(st_res, group){
  group_rank <- which(st_res$pl$data$group == group)
  gene_set <- as.character(st_res$pl$data$vertex.names[group_rank])
  enrich_res <- enrichGO(gene_set, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', pvalueCutoff = 0.1)
  enrich_pl <- dotplot(enrich_res)

  vertex.size <- 15/200
  by <- 15/200
  k <- max(st_res$plot$data$group)
  mark.col <- rainbow(length(1:k))

  pl <- st_res$plot
  pl$layers[[2]] <- NULL
  pl$layers[[2]] <- NULL
  pl <- pl + guides(size = FALSE) + guides(alpha = FALSE)
  points <- as.matrix(pl$data[which(pl$data$group==group),1:2])
  pp <- rbind(points,
              cbind(points[,1]-vertex.size-by, points[,2]),
              cbind(points[,1]+vertex.size+by, points[,2]),
              cbind(points[,1], points[,2]-vertex.size-by),
              cbind(points[,1], points[,2]+vertex.size+by))
  cl <- convex_hull(pp)
  data <- as.data.frame(cl$rescoords)
  pl$data$col <- 'Steiner genes'
  pl$data$col[group_rank] <-  'DEGs'
  pl <- pl + geom_polygon(data = data, aes(V1,V2) ,alpha = 0.3)
  pl <- pl + geom_point(aes_(x=~X1, y=~X2, color=~col), show.legend = FALSE)

  p <- plot_grid(pl, enrich_pl, nrow=2)
  res <- list(enrich_res = enrich_res, plot = p)
  return(res)
}

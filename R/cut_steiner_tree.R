#' divide the Steiner tree
#'
#' divide the Steiner tree
#'
#' @param st_res a list contains steiner tree and plot (result of PlotSteinertree)
#' @param k number of groups to divide
#'
#' @return
#' @export
#'
#' @examples
cut_steiner_tree <- function(st_res, k){
  st_tree <- st_res$tree
  comm <- cluster_edge_betweenness(st_tree)
  hc <- as.hclust(comm)
  hc_cut <- cutree(hc, k = k)
  gene_set <- list()
  for (i in 1:k) {
    gene_set[[i]] <- names(which(hc_cut == i))
  }
  pl <- st_res$plot
  pl$layers[[2]] <- NULL
  pl$layers[[2]]$mapping$colour <- NULL
  pl <- pl + guides(size = FALSE) + guides(alpha = FALSE)
  plotcord <- pl$data
  plotcord$group <- hc_cut
  plotcord$shouldLabel <- FALSE
  vertex.size <- 15/200
  by <- 15/200
  mark.col <- rainbow(length(1:k))
  data <- data.frame()
  for (i in 1:k) {
    points <- as.matrix(plotcord[which(plotcord$group==i),1:2])
    pp <- rbind(points,
                cbind(points[,1]-vertex.size-by, points[,2]),
                cbind(points[,1]+vertex.size+by, points[,2]),
                cbind(points[,1], points[,2]-vertex.size-by),
                cbind(points[,1], points[,2]+vertex.size+by))
    cl <- convex_hull(pp)
    data <- rbind(data, data.frame(cl$rescoords, group = i, col = mark.col[i]))
  }
  pl <- pl + geom_polygon(data = data, aes(x = X1, y = X2, group = group, fill = col) ,alpha = 0.3) +
    scale_fill_discrete(name = 'Group', labels=c(1:k))
  print(pl)

  st_res$plot$data <- plotcord
  return(st_res)
}

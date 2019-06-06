# Shortest Path Based Approximation (SP)
#' Title
#'
#' @param terminals
#' @param graph
#' @param color
#' @param weighted
#' @param optimize
#'
#' @return
#'
#' @examples
steinertree <- function (terminals, graph, color = T, weighted = T, optimize = T) {
  glist      <- c()
  glist[[1]] <- graph

  varlist <- check_input(terminals = terminals, glist = glist)

  glist[[1]] <- varlist[[1]]
  terminals  <- varlist[[2]]
  attr_flag  <- varlist[[3]]

  g <- glist[[1]]

  # Pick a terminal randomly and Form a subtree (sub-graph G')
  prob     <- sample(1:length(terminals), 1)
  subtree  <- terminals[[prob]]
  nsubtree <- setdiff(terminals, subtree)

  # Proceed until all terminals not in G'
  while ( !all(is.element(terminals, intersect(subtree, terminals))) ) {
    # Compute shortest paths and their lengths between each node in subtree (G') and the remaining nodes
    paths <- lapply(subtree, function (x) get.all.shortest.paths(g, x, nsubtree))

    r <- 1:length(paths)
    if(weighted == T){
      t <- sapply(r, function (r) sapply(paths[[r]][[1]], function(x)
        sum(E(g)[get.edge.ids(g,x[floor(seq(1.5,length(x),0.5))])]$weight)))
    }
    else{
      t <- sapply(r, function (r) sapply(paths[[r]][[1]], length))
    }

    # Compute a minimum for each set of lengths from each node to other nodes
    if (class(t) == "list" || class(t) == "integer") {
      r  <- 1:length(t)
      t2 <- sapply(r, function (r) min(t[[r]]))
    }
    if (class(t) == "matrix") {
      r  <- 1:dim(t)[2]
      t2 <- sapply(r, function (r) min(t[, r]))
    }

    # Find a path with minimum among minimum length
    t3 <- which(t2 == min(t2))

    # Note, graph has to have name attribute, because in found variable we assign names
    # of vertices. It is much more convenient to work with names, not with ids.
    if (length(paths) > 1) {
      if (class(t) == "list" || class(t) == "integer")
        t4 <- which(t[[t3[1]]] == min(t[[t3[1]]]))

      if (class(t) == "matrix")
        t4 <- which( t[ , t3[1]] == min(t[ , t3[1]]) )

      #found <- unlist(paths[[t3[1]]][t4][1]$res)
      found <- names(unlist(paths[[t3[1]]][t4][1]$res))
    } else {
      #found <- unlist(paths[[1]][t3][1]$res)
      found <- names(unlist(paths[[1]][t3][1]$res))
    }

    # Add all vertices from all shortest paths to subtree
    #subtree  <- union(subtree, V(g)[unique(found)])
    subtree  <- union(subtree, V(g)[unique(found)]$name)
    #nsubtree <- setdiff(nsubtree, V(g)[unique(found)])
    nsubtree <- setdiff(nsubtree, V(g)[unique(found)]$name)
  }

  # Perform "optimization": find minimum spanning tree and remove nodes of degree 1
  if (optimize) {
    steinert <- minimum.spanning.tree(induced_subgraph(g, subtree))
    a   <- V(steinert)$color
    b   <- degree(steinert, v = V(steinert), mode = c("all"))
    a1  <- match(a, "yellow")
    b1  <- match(b, "1")
    opt <- sapply(1:length(a1), function (r) a1[r] * b1[r])
    new_g <- delete.vertices(steinert, grep(1, opt))
    steinert <- new_g
  } else
    steinert <- induced_subgraph(g, subtree)

  glst <- c()
  if (color) {
    V(g)[subtree]$color   <- "green"
    V(g)[terminals]$color <- "red"

    glst[[length(glst) + 1]] <- g
  }

  glst[[length(glst) + 1]] <- steinert

  result <- glst
  if (color) {
    if (attr_flag) {
      V(result[[1]])$name <- V(result[[1]])$realname
      result[[1]] <- delete_vertex_attr(result[[1]], 'realname')
    }
  }
  if (attr_flag) {
    V(result[[length(result)]])$name <- V(result[[length(result)]])$realname
    result[[length(result)]] <- delete_vertex_attr(result[[length(result)]], 'realname')
  }

  return(result)
}

check_input <- function (terminals, glist) {

  g <- glist[[1]]
  g <- as.undirected(g)

  # Checking terminals

  if (is.null(terminals) || is.na(terminals) || length(terminals) == 0)
    stop("Error: Terminals not found")

  # Checking graph

  if (is.null(g))
    stop("Error: The graph object is Null.")

  if (length(V(g)) == 0 )
    stop("Error: The graph doesn't contain vertices.")

  if (is.null(V(g)$name)) {
    # creating name attribute
    V(g)$name <- as.character(1:length(V(g)))
    attr_flag <- FALSE
  } else {
    # creating new name and realname attributes
    V(g)$realname <- V(g)$name
    V(g)$name     <- as.character(1:length(V(g)))
    attr_flag <- TRUE
  }

  # Mathcing names of vertices and terminals, if possible

  if (class(terminals) == "character") {
    # terminals contain realname of vertices
    if (sum(terminals %in% V(g)$realname) != length(terminals)) {
      stop("Error: vertices names do not contain terminal names")
    } else {
      # Convert realnames of terminals to names (character id's)
      terminals <- V(g)$name[match(terminals, V(g)$realname)]
    }
  } else if (class(terminals) == "numeric" | class(terminals) == "integer") {
    # terminals contains id's of vertices
    terminals <- V(g)$name[terminals]
  } else
    print("Error: invalid type of terminals")

  V(g)$color            <- "yellow"
  V(g)[terminals]$color <- "red"

  varlist      <- c()
  varlist[[1]] <- g
  varlist[[2]] <- terminals
  varlist[[3]] <- attr_flag

  return(varlist)
}


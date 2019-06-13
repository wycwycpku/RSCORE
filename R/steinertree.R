#' Title get steiner tree
#'
#' @param terminals
#' @param graph
#' @param color
#' @param weighted
#' @param optimize
#' @param method
#'
#' @return
#'
#' @examples
steinertree <- function (terminals, graph, method = 'sp', color = T, weighted = T, optimize = T) {
  if(method == 'sp'){
    steinertree_sp(terminals, graph, color, weighted, optimize)
  }else if(method == 'kb'){
    steinertree_kb(terminals, graph, color, weighted, optimize)
  }else{
    stop("Provided 'method' not recognized")
  }
}

################################################
# Shortest Path Based Approximation (SP)
steinertree_sp <- function (terminals, graph, color = T, weighted = T, optimize = T) {
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

############################################################################################
# Minimum spanning tree based approximation (Kruskal's minimum spanning tree algorithm)
steinertree_kb <- function (terminals, graph, color = T, weighted = T, optimize = T) {
  makesubtrees <- function (x) {
    if ( !is.na(any(match(t3, x))) )
      #return(union(subtrees[[x]],
      #	     found[[grep(1, match(t3, x))]][[1]]))
      return(union(subtrees[[x]],
                   names(found[[grep(1, match(t3, x))]][[1]])))
    else return(subtrees[[x]])
  }
  glist      <- c()
  glist[[1]] <- graph

  varlist <- check_input(terminals = terminals, glist = glist)

  glist[[1]] <- varlist[[1]]
  terminals  <- varlist[[2]]
  attr_flag  <- varlist[[3]]

  subtreenum <- c()
  x <- c()
  g <- glist[[1]]

  # Make a Streiner Tree from every terminal
  r <- 1:length(terminals)
  subtrees  <- lapply(r, function (r) terminals[[r]])
  terminals <- subtrees
  nsubtrees <- lapply(r, function (r) setdiff(terminals, subtrees[r]))

  # Proceed until all terminals won't be added to a subtree
  while (length(subtrees) > 1) {
    # Find shortest paths between different Steiner Trees and compute their lengths
    r     <- 1:length(subtrees)
    #paths <- lapply(r, function (r) lapply(subtrees[[r]],
    #				       function (x, y) get.all.shortest.paths(g, x, y)$res,
    #				       y = nsubtrees[[r]]))
    paths <- lapply(r, function (r) lapply(subtrees[[r]],
                                           function (x, y) get.all.shortest.paths(g, x, y)$res,
                                           y = unlist(nsubtrees[[r]])))

    r <- 1:length(paths)
    if(weighted == T){
      t <- sapply(r, function (r) sapply(paths[[r]][[1]], function(x)
        sum(E(g)[get.edge.ids(g,x[floor(seq(1.5,length(x),0.5))])]$weight)))
    }
    else{
      t <- sapply(r, function (r) sapply(paths[[r]][[1]], length))
    }

    # Compute a minimum for each set of lengths from each Steiner tree to other trees
    if (class(t) == "list" | class(t) == "integer") {
      r  <- 1:length(t)
      t2 <- sapply(r, function (x) min(t[[x]]))
    }
    if (class(t) == "matrix") {
      r  <- 1:dim(t)[2]
      t2 <- sapply(r, function (r) min(t[, r]))
    }

    # Find a minimum among minimum length and paths corresponding to it
    t3    <- which(t2 == min(t2))
    t3len <- 1:length(t3)

    if (length(paths) > 1) {
      if (class(t) == "list" || class(t) == "integer" )
        t4 <- lapply(t3len, function (x) which(t[[t3[x]]] == min(t[[t3[x]]])))
      if (class(t) == "matrix")
        t4 <- lapply(t3len, function (x) which((t[ , t3[x]]) == min(t[ , t3[x]])))

      found <- lapply( t3len, function (x) paths[t3[x]][[1]][[1]][t4[[x]][1]] )
    } else {
      intersect(subtrees[[x]], V(g)[unlist(terminals)])
      print("Error")
    }

    # Merge subgraphs and paths
    subtrees <- lapply(1:length(subtrees), function (x) makesubtrees(x))

    # Delete repeated subtrees (presume that length is more than 1)
    i <- 1
    j <- 2
    while (i <= (length(subtrees) - 1)) {
      j <- i + 1
      while (j <= length(subtrees)) {
        if (length(intersect(subtrees[[i]], subtrees[[j]])) > 0) {
          subtrees[[i]] <- union(subtrees[[i]], subtrees[[j]])
          subtrees <- subtrees[-j]
          j <- j - 1
        }
        j <- j + 1
      }
      i <- i + 1
    }
    nsubtrees <- lapply(1:length(subtrees), function (x) setdiff(terminals, subtrees[[x]]))
  }

  # Perform "optimization": find minimum spanning tree and remove nodes of degree 1
  if (optimize) {
    steinert <- minimum.spanning.tree(induced_subgraph(g, subtrees[[1]]))
    a   <- V(steinert)$color
    b   <- degree(steinert, v = V(steinert), mode = c("all"))
    a1  <- match(a, "yellow")
    b1  <- match(b, "1")
    opt <- sapply(1:length(a1), function (r) a1[r] * b1[r] )
    new_g <- delete.vertices(steinert, grep(1, opt))
    steinert <- new_g
  } else
    steinert <- induced_subgraph(g, subtrees[[1]])

  glst <- c()
  if (color) {
    V(g)[subtrees[[1]]]$color     <- "green"
    V(g)[unlist(terminals)]$color <- "red"
    #V(g)[terminals]$color <- "red"

    glst[[length(glst) + 1]]  <- g
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




###############################################################
###############################################################
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


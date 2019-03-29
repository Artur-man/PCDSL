#' Drawing a triangle (2-simplex)
#'
#' @param tri_list A list of triangles (2-simplices), or a single triangle (2-simplex)
#' @param ... Additional arguements passed on to the polygon function
draw_triangle <- function(tri_list, ...) {

  if (!is.list(tri_list)){
    tri_list <- list(tri_list)
  }

  for (i in 1:length(tri_list)) {

      tri <- tri_list[[i]]

      if (!all(dim(tri)==c(3, 2))) {
        stop("The input is a not a triangle")
      }

      tri_polygon <- rbind(tri, tri[1, ])

      polygon(tri_polygon[, 1], tri_polygon[, 2], ...)

  }
}

#' Drawing a outer triangle (outer simplex in 2 dimensions)
#'
#' @param outer_tri_list A list of outer triangles, or a single outer triangle
#' @param ... Additional arguements passed on to the polygon function
draw_outer_triangle <- function(outer_tri_list,...){

  if (!is.list(outer_tri_list)) outer_tri_list <- list(outer_tri_list)

  for(i in 1:length(outer_tri_list)){

    tri <- outer_tri_list[[i]]

    if (!all(dim(tri)==c(4, 2))){
      stop("The input is a not an outer triangle")
    }

    tri_polygon <- rbind(tri[c(1,2,4,3,1), ])
    print(tri_polygon)

    polygon(tri_polygon[, 1], tri_polygon[, 2], ...)

  }
}


#' Finding the approximately minimum dominating set of graph given its adjacency matrix
#' with the greedy algorithm
#'
#' @param A The adjacency matrix
#'
#' @return S A subset (the approximate minimum dominating set) of the vertex set of the graph given the adjacency matrix A
dominate_greedy_matrix <- function(A)
{

  if (!is.matrix(A)) {
    A <- matrix(A, nrow = 1)
  }

  S <- NULL
  n <- nrow(A)
  covered <- rep(FALSE, n)

  while (!all(covered)) {
    od <- apply(A, 1, sum)
    i <- which.max(od)
    covered[A[i, ]==TRUE] <- TRUE
    S <- c(S, i)
    A[, covered==TRUE] <- FALSE
  }

  return(S)
}

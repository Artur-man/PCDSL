#' Drawing a triangle
#'
#' @param tri_list A list of triangles, or a single triangle
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

#' Drawing an outer triangle
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
#' @return S A subset of the vertex set of the graph given the adjacency matrix A
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

#' Finding the exact minimum dominating set of graph given its adjacency matrix
#' with the brute force algorithm (checking all possible solutions)
#'
#' @param A The adjacency matrix
#'
#' @return S A subset of the vertex set of the graph given the adjacency matrix A
dominate_brute_matrix <- function(A)
{

  if (!is.matrix(A)) {
    A <- matrix(A, nrow=1)
  }

  n <- nrow(A)
  D <- NULL
  for (i in 1:n) {
    cmb <- combn(n, i)
    for( k in 1:ncol(cmb)) {
      sub_A <- A[cmb[, k], ]
      sub_A <- matrix(sub_A, ncol = n)
      mat_check <- apply(sub_A, 2, any)
      if (all(mat_check)) {
        D <- cmb[, k]
        break
      }
    }
    if (!is.null(D)) {
      break
    }
  }
  return(D)
}

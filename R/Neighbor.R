#' Generate Random Neighboring Points
#'
#' This function generates a vector of random neighboring points around a given point `x`, constrained by the range specified in matrix `D`.
#'
#' @param x A numeric vector representing the initial point.
#' @param D A numeric matrix with two columns. Each row represents the minimum and maximum bounds for the corresponding element of `x`.
#'
#' @return A numeric vector of the same length as `x`, containing the generated neighboring points.
#' @export
#'
#' @examples
#' x <- c(1, 2)
#' D <- matrix(c(0, 3, 1, 4), ncol = 2, byrow = TRUE)
#' Neighbor(x, D)
Neighbor <- function(x, D) {
  x <- matrix(x, nrow = 1)
  m <- x - (x - D[, 1]) / 3
  M <- x + (D[, 2] - x) / 3
  xv <- stats::runif(ncol(x), min = m, max = M)
  return(xv)
}


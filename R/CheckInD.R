#' @title Constraint Handling Function for Optimization
#'
#' @description
#' The `TT` function ensures that all values in a matrix `X` are constrained
#' within specified bounds defined by a matrix `D`.
#'
#' @param X A numeric matrix where each row represents a candidate solution in the optimization process.
#' @param D A numeric matrix with two columns wich represent the domain.
#'
#' @return A matrix `X` with the same dimensions as the input, but with all values
#' constrained within the bounds defined by `D`.
#'
#' @details
#' This function iterates through each element of the matrix `X`. For each element,
#' it checks if the value is within the bounds specified by the corresponding row in `D`.
#' If the value is below the lower bound, it is replaced by the lower bound.
#' If it is above the upper bound, it is replaced by the upper bound.
#' @examples
#' X <- matrix(runif(10), nrow=5, ncol=2)
#' D <- matrix(c(0.2, 0.8, 0.1, 0.9), nrow=2, ncol=2, byrow=TRUE)  # Bounds
#' Xc <- TT(X, D)
#' print(Xc)
#'
#' @export
TT <- function(X, D) {
  N <- nrow(X)
  for(i in 1:N) {
    for(j in 1:nrow(D)) {
      X[i,j] <- min(max(X[i,j], D[j,1]), D[j,2])
    }
  }
  return(X)
}

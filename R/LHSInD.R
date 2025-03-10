#' Latin Hypercube Sampling in domain D
#'
#' Generates a Latin Hypercube Sample (LHS) scaled to fit within the domain.
#'
#' @param N The number of samples to generate.
#' @param D A matrix where each row represents the range min, max for a design variable.
#' @return A matrix of size N x d, where d is the number of design variables (rows of D).
#' @import lhs
#' @export
#' @examples
#' D <- matrix(c(0, 10, 0, 5), ncol = 2, byrow = TRUE)
#' LHSD(100, D)
LHSD <- function(N,D){
  d <- nrow(D)
  G <- lhs::maximinLHS(N,d)

  for(i in 1:d){
    G[,i] <- (D[i,2]-D[i,1])*G[,i] + D[i,1]
  }

  G = as.matrix(G)

  return(G)
}

#' Gridoptim: Grid-Based Optimization Function
#'
#' This function performs an optimization using a grid-based search approach. It iteratively selects the best points from a Latin Hypercube Sampling (LHS) design, optimizing the objective function provided.
#'
#' @param fobj A function representing the objective function to be optimized. It should accept a matrix as input where each row is a different set of parameters.
#' @param size An integer indicating the number of points to select during the optimization process.
#' @param D A matrix with two columns, specifying the lower and upper bounds for each dimension of the design space. Each row corresponds to one dimension.
#' @param nLHS An integer specifying the number of samples to generate in the Latin Hypercube Sampling (default is 1000).
#' @param isMax A logical value indicating whether to maximize (`TRUE`, default) or minimize (`FALSE`) the objective function.
#'
#' @return A list with two elements:
#' \item{Xopt}{A matrix of the selected optimal points. Each row represents a point in the design space.}
#' \item{Copt}{The optimal value of the objective function at the selected points.}
#'
#' @details The function starts by generating a large number of points using Latin Hypercube Sampling (LHS). It then iteratively selects the best point according to the objective function, and removes that point from the set of candidates. This process is repeated until the specified number of points (`size`) is selected.
#'
#' @examples
#' \dontrun{
#' # Define an example objective function
#' fobj <- function(x) sum(x^2)
#'
#' # Define the design space
#' D <- matrix(c(-5, 5, -5, 5), ncol=2, byrow=TRUE)
#'
#' # Run the optimization
#' R <- Gridoptim(fob = fobj, size = 10, D = D, nLHS = 1000, isMax = FALSE)
#'
#' # View the results
#' print(R$Xopt)
#' print(R$Copt)
#' }
#' @export
Gridoptim <- function(fobj, size, D, nLHS=1000, isMax=TRUE){

  if (!is.numeric(size)) stop("size must be numeric")
  if (!is.matrix(D) || ncol(D) != 2) stop("D must be a 2-column matrix")
  if (!is.numeric(nLHS)) stop("nLHS must be numeric")
  if (!is.logical(isMax)) stop("isMax must be logical (TRUE or FALSE)")

  XX = LHSD(nLHS, D)
  Xopt = matrix(0, ncol=nrow(D), nrow=size)
  Xopt[1,] = stats::runif(nrow(D), min=D[, 1], max=D[, 2])
  for(i in 2:size){
    fXX = apply(XX, 1, function(x){return(fobj(rbind(Xopt[1:(i-1),], x)))})

    if(isMax){istar = which.max(fXX)}else{istar = which.min(fXX)}

    Xopt[i, ] <- XX[istar, ]
    XX = matrix(XX[-istar, ], ncol=nrow(D))
    if(i==size) Copt = fXX[istar]
  }

  R = list(Xopt = Xopt, Copt=Copt)
  return(R)
}

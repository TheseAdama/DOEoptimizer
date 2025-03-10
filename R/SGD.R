#' @title Stochastic Gradient Descent Optimization for Design of Experiments
#' @description Implements a stochastic gradient descent (SGD) optimization algorithm for
#' optimizing a function in the context of design of experiments (DOE). The algorithm iteratively
#' updates the design to either maximize or minimize the objective function.
#'
#' @param fobj A function representing the objective function to be optimized. It should accept a matrix as input where each row is a different set of design.
#' @param D A numeric matrix with two columns defining the bounds of the design to be optimized.
#' @param size An integer specifying the size of the experimental design (number of points).
#' @param maxiter An integer specifying the maximum number of iterations (default: 10000).
#' @param Xinit An optional numeric matrix specifying the initial design. If \code{NULL}, a design will be generated using Latin Hypercube Sampling (LHS).
#' @param a A numeric value for the initial step size (default: 1e-3).
#' @param c A numeric value for the coefficient controlling the amplitude of the perturbations (default: 0.01).
#' @param A An optional numeric constant used for adjusting the step size in the geometric decay (default: \code{NULL}, computed as 0.1 * \code{maxiter}).
#' @param gamma A numeric value between 0 and 1 for the exponent controlling the decay rate of the perturbation amplitude (default: 1/6).
#' @param isMax A logical value indicating whether the objective function should be maximized (\code{TRUE}) or minimized (\code{FALSE}).
#'
#' @details The algorithm uses a stochastic approach to estimate the gradient of the objective function and iteratively updates the design in the direction of the gradient. The step size diminishes over time, following a geometric decay to ensure convergence.
#' This approach is based on the Simultaneous Perturbation Stochastic Approximation (SPSA) method, which is particularly efficient for high-dimensional optimization problems. The method is well-documented in the following references:
#'
#' - Huan, Xun and Marzouk, Youssef M. (2013). Simulation-based optimal Bayesian experimental design for nonlinear systems. *Journal of Computational Physics*, 232, 288-317.
#' - Spall, J. C. (1998). An overview of the simultaneous perturbation method for efficient optimization. *Johns Hopkins APL Technical Digest*, 19, 482-492.
#' - Kleinman, N. L., Spall, J. C., & Naiman, D. Q. (1999). Simulation-based optimization with stochastic approximation using common random numbers. *Management Science*, 45, 1570-1578.
#' - Spall, J. C. (1992). Multivariate stochastic approximation using a simultaneous perturbation gradient approximation. *IEEE Transactions on Automatic Control*, 37, 332-341.
#'
#' @return A list with optimized design \code{Xopt} and the corresponding objective value \code{Copt}.
#' @examples
#' # Example usage of SGDoptim
#' fobj <- function(X){return(log(det(t(X)%*%X)))}
#' D <- matrix(c(0, 10, 0, 5), ncol = 2, byrow = TRUE)
#' R <- SGDoptim(fobj, D, size=10, maxiter=1000)
#' print(R$Xopt)
#' print(R$Copt)
#'
#' @export
SGDoptim <- function(fobj, D, size, maxiter=NULL, Xinit=NULL,
                     a=1e-4 , c=0.01, A=NULL, gamma=1/6, isMax = TRUE) {

  if(is.null(maxiter)){
    maxiter <- if(nrow(D) <= 4) {10^(2 + nrow(D)) } else { 10^6 }
  }
  if (!is.numeric(maxiter)) stop("maxiter must be numeric")
  if (!is.numeric(size)) stop("size must be numeric")
  if (!is.numeric(a)) stop("a must be numeric")
  if (!is.numeric(c) || (c > 1 || c < 0)) stop("c must be a numeric between 0 and 1")
  if (!is.numeric(gamma) || (gamma > 1 || gamma < 0)) stop("gamma must be numeric between 0 and 1")
  if (!is.logical(isMax)) stop("isMax must be logical (TRUE or FALSE)")
  if (!is.matrix(D) || ncol(D) != 2) stop("D must be a 2-column matrix")
  if (!is.null(Xinit) && !is.matrix(Xinit)) stop("Xinit must be a matrix if provided")

  if(is.null(A)) { A <- 0.1 * maxiter }
  if (!is.numeric(A)) stop("A must be numeric")
  c <- c * min(D[, 2] - D[, 1])

  if(is.null(Xinit)) {
    X <- LHSD(size, D)
  } else {
    X <- Xinit
  }
  hh = rep(0,maxiter)
  hh[1] = fobj(X)
  for (k in 2:maxiter) {
    DD <- matrix(2 * stats::rbinom(nrow(D) * size, 1, 0.5) - 1, ncol=nrow(D))
    ck <- c / (k^gamma)

    Xa <- X + ck * DD
    Xb <- X - ck * DD
    Gk <- ((fobj(Xa) - fobj(Xb)) / (2 * ck)) * DD

    ak <- a*(1-1e-6)**k #a / (k + A)^gamma

    if(isMax) {
      X <- X + ak * Gk
    } else {
      X <- X - ak * Gk
    }

    X <- TT(X, D)

    hh[k] = fobj(X)
    if (k %% 100 == 0) {
      pp <- (k / maxiter) * 100
      cat(sprintf("Progress: %.2f%%\n", pp))
      Sys.sleep(0.001)
    }
  }
  plot(c(1:maxiter), hh, main="History SGD",xlab="Iteration",
       ylab="fobj", lwd=1.5, type='l', col="blue",cex.axis=1.5, cex.lab=1.5, cex.main=1.8)

  R = list(Xopt = X, Copt=fobj(X))
  return(R)
}

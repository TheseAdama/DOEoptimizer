#' FIBoptim: Forward, Improvement, and Backward Optimization Method
#'
#' This function implements a three-step optimization process-Forward, Improvement, and Backward (FIB)-to optimize a design criterion.
#' The method is applicable provided that the function is either increasing (XcY => fobj(X)<=fobj(X)) for maximization problem
#' or decreasing (X c Y => fobj(X)<=fobj(X)) for minimization problem.
#' @param fobj A function representing the objective function to be optimized.
#' @param D A numeric matrix where each row represents the range min, max for a design variable.
#' @param size The target size for the final design (number of rows).
#' @param nF The number of forward iterations. Default is 100.
#' @param nI The number of improvement iterations. Default is 100.
#' @param nLHS The number of Latin Hypercube Samples (LHS) to generate. Default is 1000.
#' @param isMax Logical. If TRUE, the algorithm maximizes the objective function; if FALSE, it minimizes it. Default is TRUE.
#' @return An updated DOEoptim object with the optimal design (`Xopt`) and corresponding criterion value (`Copt`).
#' @export
#'
#' @examples
#' fobj <- function(X){return(det(t(X)%*%X))}
#' D <- matrix(c(0, 10, 0, 5), ncol = 2, byrow = TRUE)
#' R <- FIBoptim(fobj, D, size=10, nF=100, nI=100, nLHS=500, isMax = TRUE)
#'
#' print(R$Xopt)
#' print(R$Copt)

FIBoptim <- function(fobj, D, size, nF = 100, nI = 100, nLHS = 1000, isMax = TRUE) {

  if (!is.matrix(D) || ncol(D) != 2) stop("D must be a 2-column matrix")
  if (!is.numeric(size) || length(size) != 1) stop("size must be a single numeric value")
  if (!is.numeric(nF) || length(nF) != 1) stop("nF must be a single numeric value")
  if (!is.numeric(nI) || length(nI) != 1) stop("nI must be a single numeric value")
  if (!is.numeric(nLHS) || length(nLHS) != 1) stop("nLHS must be a single numeric value")
  if (!is.logical(isMax) || length(isMax) != 1) stop("isMax must be a single logical value")

  if (nF > nLHS) warning("nF should be less than or equal to nLHS")
  XX <- LHSD(nLHS, D)
  XF <- LHSD(1, D)
  hh = rep(0, nF)
  hh[1] = fobj(XF)
  for (i in 2:nF) {
    cXX <- apply(XX, 1, function(x) {return(fobj(rbind(XF, x)))})
    istar <- if (isMax) which.max(cXX) else which.min(cXX)
    XF <- rbind(XF, XX[istar, ])
    XX <- XX[-istar, ]
    XX <- matrix(XX, ncol = nrow(D))
    hh[i] = fobj(XF)

    if (i %% 10 == 0) {
      cat(sprintf("Forward step: %.2f%%\n", (i / nF) * 100))
      Sys.sleep(0.001)
    }

  }
  plot(c(1:nF), hh, main="History Forward Step",xlab="Iteration",
       ylab="fobj", lwd=1.5, type='l', col="blue", cex.axis=1.5, cex.lab=1.5, cex.main=1.8)
  Cref <- fobj(XF)
  XX <- LHSD(nLHS, D)
  hh <- rep(0, nI)
  for (j in 1:nI) {
    ii <- sample(1:nF, 1)
    cXX <- apply(XX, 1, function(x){
      Xc <- XF
      Xc[ii, ] <- x
      return(fobj(Xc))
    })
    istar <- if (isMax) which.max(cXX) else which.min(cXX)
    Ccurrent <- cXX[istar]
    if (Cref <= Ccurrent) {
      XF[ii, ] <- XX[istar, ]
      XX <- XX[-istar, ]
      XX <- matrix(XX, ncol = nrow(D))
      Cref = Ccurrent
    }
    hh[j] = Cref
    if (j %% 10 == 0) {
      cat(sprintf("Improvement step: %.2f%%\n", (j / nI) * 100))
      Sys.sleep(0.001)
    }

  }
  plot(c(1:nI), hh, main="History Improvement Step",xlab="Iteration",
       ylab="fobj", lwd=1.5, type='l', col="green",cex.axis=1.5, cex.lab=1.5, cex.main=1.8)

  hh <- rep(0, nF - size)
  for (k in 1:(nF - size)) {
    cXF <- apply(array(1:nrow(XF)), 1, function(ii) {return(fobj(matrix(XF[-ii, ], ncol = nrow(D))) ) } )
    istar <- if (isMax) which.max(cXF) else which.min(cXF)
    XF <- matrix(XF[-istar, ], ncol = nrow(D))
    if (k %% 10 == 0) {
      cat(sprintf("Backward step: %.2f%%\n", (k / (nF - size)) * 100))
      Sys.sleep(0.001)
    }
    hh[k] = fobj(XF)
  }
  plot(c(1:(nF-size)), hh, main="History Backward Step",xlab="Iteration",
       ylab="fobj", lwd=1.5, type='l', col="red",cex.axis=1.5, cex.lab=1.5, cex.main=1.8)

  R = list(Xopt = XF, Copt=fobj(XF))
  return(R)
}


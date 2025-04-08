#' Simulated Annealing Method for Optimal Design of Experiments
#'
#' This function optimizes a criterion using the Simulated Annealing (SA) method. It iteratively searches for the best experimental design within a specified domain by probabilistically exploring neighboring designs.
#'
#' @param fobj A function representing the objective function to be optimized. It should accept a matrix as input where each row is a different set of experiments.
#' @param size An integer specifying the size of the experimental design (number of points).
#' @param D A numeric matrix representing the domain of design elements. Each row defines the range for one design variable.
#' @param maxiter Maximum number of iterations for the algorithm (default: 10^(1 + nrow(D)).
#' @param temp Initial temperature for the annealing process (default: 0.1).
#' @param c Cooling rate, controlling the reduction of temperature over iterations (default: 0.99).
#' @param Xinit Initial design matrix for the experiments. If NULL, the initial design from `object` is used.
#' @param ptype The probability scheme used for accepting new designs. Can be "metropolis" (default) or "glauber".
#' @param tempscheme The temperature scheme used to reduce the temperature. Can be "geo" for geometric (default) or "log" for logarithmic.
#' @param isMax The parameter that specifies whether the function is to be maximized or minimized. Its default value is "TRUE".
#' @return A list with the optimal design (`Xopt`) and the corresponding criterion value (`Copt`) found.
#' @export
#'
#' @examples
#' fobj <- function(X){return(log(det(t(X)%*%X)))}
#' D <- matrix(c(0, 10, 0, 5), ncol = 2, byrow = TRUE)
#' R <- SAoptim(fobj, size=10, D=D, maxiter = 1000, temp = 1, c = 0.95)
#' print(R$Xopt)
#' print(R$Copt)
SAoptim<- function(fobj, size, D, maxiter = NULL, temp = 0.1, c = 0.99, Xinit = NULL,
                   ptype = "metropolis", tempscheme = "geo", isMax = TRUE) {

  if(is.null(maxiter)){maxiter <- if(nrow(D) <= 4) {10^(1 + nrow(D))} else { 10^6 }}
  if (is.null(Xinit)){Xinit = LHSD(N=size,D)}
  if (!is.matrix(D)) stop("D must be a 2-column matrix")
  if (!is.matrix(Xinit)) stop("Xinit must be a matrix")
  if (ncol(Xinit) != nrow(D)) stop("Xinit must have the same number of columns as the number of rows in D")
  if (!is.numeric(maxiter)) stop("maxiter must be numeric")
  if (!is.numeric(temp)) stop("temp must be numeric")
  if (!is.character(ptype)) stop("ptype must be a character string")
  if (!is.character(tempscheme)) stop("tempscheme must be a character string")

  if (!(ptype %in% c("metropolis", "glauber"))) stop("Invalid ptype. Use 'metropolis' or 'glauber'.")
  if (!(tempscheme %in% c("geo", "log"))) stop("Invalid tempscheme. Use 'geo' or 'log'.")
  if(!is.logical(isMax)) stop("IsMax must be logical (TRUE or FALSE)")

  Xopt <- Xinit
  Copt <- fobj(Xopt)
  hh <- rep(0, maxiter)
  hh[1] = Copt
  for (k in 2:maxiter) {
    for (i in 1:nrow(Xopt)) {
      Xtemp <- Xopt
      Xtemp[i, ] <- Neighbor(Xopt[i, ], D)  # Generate neighboring design
      Ck <- fobj(Xtemp)

      if (is.na(Ck)) next

      if(isMax){dd <- Copt - Ck}else{dd <- Ck - Copt}

      if (ptype == "metropolis") {aa <- min(1, exp(-dd / temp))
      } else if (ptype == "glauber") {aa <- exp(-dd / temp) / (1 + exp(-dd / temp))}

      if (stats::runif(1) < aa) {
        Xopt <- Xtemp
        Copt <- Ck}
    }
    # Update temperature based on the chosen scheme
    if (tempscheme == "geo") {
      temp <- temp * c
    } else if (tempscheme == "log") {
      temp <- temp / log(1 + k)
    }

    # Affichage de la progression
    if (k %% 100 == 0) {
      pp <- (k / maxiter) * 100
      cat(sprintf("Progress: %.2f%%\n", pp))
      Sys.sleep(0.001)
    }
    hh[k] = fobj(Xopt)
  }
  plot(c(1:maxiter), hh, main="History SA",xlab="Iteration",
       ylab="fobj", lwd=1.25, type='l', col="blue", cex.axis=1.5, cex.lab=1.5, cex.main=1.8)

  R = list(Xopt = Xopt, Copt = Copt)
  return(R)
}

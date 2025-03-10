#' DOEoptimizer: A Package for Optimizing Experimental Design Criteria
#'
#' \code{DOEoptimizer} is a package that implements four optimization algorithms
#' specifically designed for optimizing experimental design criteria. The package provides tools to
#' efficiently generate and optimize design matrices according to various criteria commonly used in
#' design of experiments (DOE).
#'
#' @details
#' The \code{DOEoptimizer} package includes the following optimization algorithms:
#' \itemize{
#'   \item \strong{GEAoptim}: Genetic Algorithm Optimization for DOE criteria.
#'   \item \strong{SAoptim}: Simulated Annealing Optimization for DOE criteria.
#'   \item \strong{FIBoptim}: Forward, Improvement, and Backward greedy Optimization Algorithm for DOE criteria.
#'   \item \strong{SGDoptim}: Stochastic Gradient Descent Optimization for DOE criteria.
#'   \item \strong{Gridoptim}: Grid-Based algorithm
#' }
#'
#' The package leverages existing R packages like \code{stats}, \code{base}, and \code{lhs} to perform
#' various calculations and sampling tasks. The main objective is to provide a flexible and powerful
#' set of tools for researchers and engineers working in the field of experimental design.
#'
#' @note
#' This package is released under the \code{GPL-3} license.
#'
#' @references
#' Barry, A. (2024). \emph{DOEoptimizer: Algorithms for Experimental Design Optimization}. R Package Version 1.0.
#'
#' @docType package
#' @name DOEoptimizer-package
"_PACKAGE"

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Grid.R
\name{Gridoptim}
\alias{Gridoptim}
\title{Gridoptim: Grid-Based Optimization Function}
\usage{
Gridoptim(fobj, size, D, nLHS = 1000, isMax = TRUE)
}
\arguments{
\item{fobj}{A function representing the objective function to be optimized. It should accept a matrix as input where each row is a different set of parameters.}

\item{size}{An integer indicating the number of points to select during the optimization process.}

\item{D}{A matrix with two columns, specifying the lower and upper bounds for each dimension of the design space. Each row corresponds to one dimension.}

\item{nLHS}{An integer specifying the number of samples to generate in the Latin Hypercube Sampling (default is 1000).}

\item{isMax}{A logical value indicating whether to maximize (\code{TRUE}, default) or minimize (\code{FALSE}) the objective function.}
}
\value{
A list with two elements:
\item{Xopt}{A matrix of the selected optimal points. Each row represents a point in the design space.}
\item{Copt}{The optimal value of the objective function at the selected points.}
}
\description{
This function performs an optimization using a grid-based search approach. It iteratively selects the best points from a Latin Hypercube Sampling (LHS) design, optimizing the objective function provided.
}
\details{
The function starts by generating a large number of points using Latin Hypercube Sampling (LHS). It then iteratively selects the best point according to the objective function, and removes that point from the set of candidates. This process is repeated until the specified number of points (\code{size}) is selected.
}
\examples{
\dontrun{
# Define an example objective function
fobj <- function(x) sum(x^2)

# Define the design space
D <- matrix(c(-5, 5, -5, 5), ncol=2, byrow=TRUE)

# Run the optimization
R <- Gridoptim(fob = fobj, size = 10, D = D, nLHS = 1000, isMax = FALSE)

# View the results
print(R$Xopt)
print(R$Copt)
}
}

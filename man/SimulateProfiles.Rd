% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimulateProfiles.R
\name{SimulateProfiles}
\alias{SimulateProfiles}
\title{Simulate profiles}
\usage{
SimulateProfiles(phytree, alpha, B, root=F)
}
\arguments{
\item{phytree}{a phylogenetic tree.}

\item{alpha}{a vector of the intrinsic rates (\eqn{\alpha}) in CCM. The length of the vector decides how many profiles to simulate.}

\item{B}{a symmetric association matrix (\eqn{B}) in CCM. If the input matrix is not symmetric, the upper triangle of the matrix will be used.}

\item{root}{a vector of states at root node. `root=F` as default (random states will be assigned at root).}
}
\value{
a matrix containing simulated profiles.
}
\description{
Simulate binary profiles based on the Community Coevolution Model
using the input tree and user-defined parameters.
}
\examples{
set.seed(123)
# generate a random 50-tip tree
t <- rtree(50)
# simulate 5 profiles
n <- 5
# assigning parameter values.
# The parameters should be in a reasonable scale otherwise the simulated profiles may have all 0s or 1s.
alpha <- runif(n, -0.5, 1) # assign n random intrinsic rates.
association matrix B
B <- matrix(0,n,n)
B[1:3,1:3] <- 1 # first 3 genes are correlated. The other 2 have no interaction with others
diag(B) <- runif(n, -0.5,0.5) # assign half the difference between gain and loss rates for each gene
# simulate 5 random profiles
simDF <- SimulateProfiles(t, alpha, B)
# plot the profiles
d <- TreeToDend(t) # convert the tree to dendrogram
ProfilePlot(simDF, d)
}

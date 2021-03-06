% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.functions.R
\name{permiv}
\alias{permiv}
\title{Permutation of group membership for a single candidate partition}
\usage{
permiv(grp, spp, ivTot = ivTot, nPerm = 250)
}
\arguments{
\item{grp}{A vector of binary (0,1) assignments indicating group membership
for a single partition of nrows in 'spp' according to a ranks of an
environmental variable.}

\item{spp}{A site-by-taxon matrix of taxa counts at each sampling location.}

\item{ivTot}{A logical indicating whether IndVal scores should be calculated
using total relative abundance or the mean relative abundace originally
proposed by Dufrene and Legendre (1997).}

\item{nPerm}{The number of permutations to be performed.}
}
\value{
A matrix (ivmax) of IndVal maxima with nrow equal to the number of
permutations and ncol equal to the number of taxa in 'spp'.
}
\description{
Randomizes group assignments for all permutations based on ranked values of
the environmental gradient, then calls the appropriate matrix-based IndVal
function.
}
\details{
This function handles the randomization portion of the permutation procedure,
and then serves as a wrapper for 'indvalp' and 'indvalps' (depending on the
value of 'ivTot'), which estimate change-point distributions across all
permutations.  The output is the same as 'indvalp' and 'indvalps' because the
function simply passes their products on.
}
\references{
Dufrene, M. and P. Legendre. 1997. Species assemblages and
indicator species: the need for a flexible asymmetrical approach. Ecol.
Mon. 67:345-366.

Baker, ME and RS King.  2010. A new method for detecting and
interpreting biodiversity and ecological community thresholds. Methods in
Ecology and Evolution 1(1): 25:37.
}
\seealso{
\code{\link[=indvalp]{indvalp()}}, \code{\link[=indvalps]{indvalps()}}, \code{\link[=getivz]{getivz()}}
}
\author{
M. Baker and R. King
}
\keyword{TITAN}
\keyword{permutation}

\name{calculate.and.compare.correlations}
\alias{calculate.and.compare.correlations}
\title{Calculate and compare two correlations}
\description{Calculate two correlation coefficients and test if they are statistically different (given the sample-sizes).}
\usage{
calculate.and.compare.correlations(x1, x2, y1, y2, method = 'spearman');
}
\arguments{
  \item{x1}{Vector of x-values for correlation 1}
  \item{x2}{Vector of x-values for correlation 2}
  \item{y1}{Vector of y-values for correlation 1}
  \item{y2}{Vector of y-values for correlation 2}
  \item{method}{What correlation method should be used? Currently only supports Spearman's rho}
  }
\value{No return value}
\author{Paul C. Boutros}
\note{I worry about the use of length(x1) and length(x2) in the code. We really need to be determining the number of jointly non-NA entries (since we set use = 'pairwise.complete.obs'). I think this would be something like length( x1[ !is.na(x1&y1) ] );}
\examples{
calculate.and.compare.correlations(
	x1 = rnorm(100),
	x2 = rnorm(100),
	y1 = runif(100),
	y2 = runif(100),
	method = 'spearman'
	);
}
\seealso{
  \code{\link{compare.correlations}}
  }
\keyword{htest}

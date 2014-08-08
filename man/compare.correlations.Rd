\name{compare.correlations}
\alias{compare.correlations}
\title{Compares two Correlations}
\description{Tests two correlation coefficients (given the sample-sizes) to determine if they are statistically similar.}
\usage{
compare.correlations(corr1, corr2, n1, n2, method);
}
\arguments{
  \item{corr1}{First Spearman's Rho (raw)}
  \item{corr2}{Second Spearman's Rho (raw)}
  \item{n1}{First sample-size}
  \item{n2}{Second sample-size}
  \item{method}{What correlation method should be used? Currently only supports Spearman's rho}
}
\value{Returns a vector of length nine containing: corr1, corr2, sd1, sd2, delta.corr, pooled SD, p-value, n1, n2}
\references{For Spearman's rho we use the formulae from: Fieller et al. "Tests for Rank Correlation Coefficients. I," Biometrica, Vol. 44, No. 3/4 (Dec. 1957)}
\author{Mehrdad Shamsi}
\examples{
# this small differential correlation is actually borderline significant (p = 0.10) 
compare.correlations(
	corr1 = 0.15,
	corr2 = -0.15,
	n1 = 55,
	n2 = 80,
	method = 'spearman'
	);
}
\seealso{
  \code{\link{calculate.and.compare.correlations}}
  }
\keyword{htest}

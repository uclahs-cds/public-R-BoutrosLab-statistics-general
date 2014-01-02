\name{get.utest.p.and.medianfoldchange}
\alias{get.utest.p.and.medianfoldchange}
\title{Get the p-value from a u-test and calculate the median foldchange}
\description{Tests the median difference of two samples and returns the p-value from a u-test with useful default parameters and calculates the median foldchange. Returns NA if the u-test is uncalculable (i.e. does try() handling).}
\usage{
get.utest.p.and.medianfoldchange(x, group1, group2, paired = FALSE, alternative = 'two.sided', logged = TRUE);
}
\arguments{
  \item{x}{Vector of numbers to analyze}
  \item{group1}{Vector of TRUE/FALSE indicating which are the x samples for the u-test}
  \item{group2}{Vector of TRUE/FALSE indicating which are the y samples for the u-test}
  \item{paired}{Do a paired u-test?}
  \item{alternative}{What's the null-hypothesis?}
  \item{logged}{Indicate whether data is log-transformed or not}
}
\details{Calculates fold-changes as median(group2)/median(group1) in normal-space and median(group2) - median(group1) in log-space, with appropriate NA-handling.}
\value{Returns a two-element vector containing the p-value and fold-change, or NAs if an error occurs.}
\author{Cindy Q. Yao}
\examples{
get.utest.p.and.medianfoldchange(
	x = rnorm(100),
	group1 = c( rep(TRUE, 50), rep(FALSE, 50) ),
	group2 = c( rep(FALSE, 50), rep(TRUE, 50) ),
	logged = FALSE
	);
}
\seealso{
  \code{\link{get.utest.p}}
  }
\keyword{htest}

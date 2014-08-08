\name{get.ttest.p.and.foldchange}
\alias{get.ttest.p.and.foldchange}
\title{Get the p-value from a t-test and calculate the foldchange}
\description{Returns the p-value from a t-test with useful default parameters and calculates the foldchange. Returns NA if the t-test is uncalculable (i.e. does try() handling).}
\usage{
get.ttest.p.and.foldchange(x, group1, group2, paired = FALSE, var.equal = FALSE, alternative = 'two.sided', logged = TRUE);
}
\arguments{
  \item{x}{Vector of numbers to analyze}
  \item{group1}{Vector of TRUE/FALSE indicating which are the x samples for the t-test}
  \item{group2}{Vector of TRUE/FALSE indicating which are the y samples for the t-test}
  \item{paired}{Do a paired t-test?}
  \item{var.equal}{Assume Welch's correction?}
  \item{alternative}{What's the null-hypothesis?}
  \item{logged}{Indicate whether data is log-transformed or not}
}
\details{Calculates fold-changes as mean(group2)/mean(group1) in normal-space and mean(group2) - mean(group1) in log-space, with appropriate NA-handling.}
\value{Returns a two-element vector containing the p-value and fold-change, or NAs if an error occurs.}
\author{Paul C. Boutros}
\examples{
get.ttest.p.and.foldchange(
	x = rnorm(100),
	group1 = c( rep(TRUE, 50), rep(FALSE, 50) ),
	group2 = c( rep(FALSE, 50), rep(TRUE, 50) ),
	logged = FALSE
	);
}
\seealso{
  \code{\link{get.ttest.p}}
  }
\keyword{htest}

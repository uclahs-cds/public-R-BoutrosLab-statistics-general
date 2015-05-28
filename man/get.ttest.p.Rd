\name{get.ttest.p}
\alias{get.ttest.p}
\title{Get the p-value from a t-test}
\description{Returns the p-value from a t-test with useful default parameters. Returns NA if the t-test is uncalculable (i.e. does try() handling).}
\usage{
get.ttest.p(x, group1, group2, paired = FALSE, var.equal = FALSE, 
    alternative = 'two.sided');
}
\arguments{
  \item{x}{Vector of numbers to analyze}
  \item{group1}{Vector of TRUE/FALSE indicating which are the x samples for the t-test}
  \item{group2}{Vector of TRUE/FALSE indicating which are the y samples for the t-test}
  \item{paired}{Do a paired t-test?}
  \item{var.equal}{Assume Welch's correction?}
  \item{alternative}{What's the null-hypothesis?}
}
\value{No return value}
\author{Paul C. Boutros}
\examples{
get.ttest.p(
	x = rnorm(100),
	group1 = c(rep(TRUE, 50), rep(FALSE, 50)),
	group2 = c(rep(FALSE, 50), rep(TRUE, 50)),
	paired = FALSE,
	var.equal = FALSE,
	alternative = 'two.sided'
	);
}
\seealso{
  \code{\link{get.ttest.p.and.foldchange}}
  }
\keyword{htest}

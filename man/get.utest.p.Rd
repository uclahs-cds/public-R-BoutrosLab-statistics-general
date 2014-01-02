\name{get.utest.p}
\alias{get.utest.p}
\title{Get the p-value from a u-test}
\description{Performs a non-parametric test by assigning ranks to individual observations in the samples. Tests the median differences between the compared samples. Returns the p-value from a u-test with useful default parameters. Returns NA if the u-test is uncalculable (i.e. does try() handling).}
\usage{
get.utest.p(x, group1, group2, paired = FALSE, alternative = 'two.sided');
}
\arguments{
  \item{x}{Vector of numbers to analyze}
  \item{group1}{Vector of TRUE/FALSE indicating which are the x samples for the u-test}
  \item{group2}{Vector of TRUE/FALSE indicating which are the y samples for the u-test}
  \item{paired}{Do a paired u-test?}
  \item{alternative}{What's the null-hypothesis?}
}
\details{If group2 = NULL, or if both group1 and group2 are given and paired is TRUE, the Wilcoxon signed rank test will be performed. If both group1 and group2 are given, the Wilcoxon rank sum test will be performed.}
\value{returns the p value from u-test between the mean of group1 and group2}
\author{Cindy Q. Yao}
\examples{
get.utest.p(
	x = rnorm(100),
	group1 = c(rep(TRUE, 50), rep(FALSE, 50)),
	group2 = c(rep(FALSE, 50), rep(TRUE, 50)),
	paired = FALSE,
	alternative = 'two.sided'
	);
}
\seealso{
  \code{\link{get.utest.p.and.medianfoldchange}}
  }
\keyword{htest}

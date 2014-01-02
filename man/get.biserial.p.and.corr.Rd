\name{get.biserial.p.and.corr}
\alias{get.biserial.p.and.corr}
\title{Calculate a point-biserial correlation and its statistical significance}
\description{Returns the point-biserial correlation and p-value for two variables (one being dichotomous) using a user-specified correlation metric. P-values are estimated analytically, not via permutation-testing.}
\usage{
get.biserial.p.and.corr(x, y, alternative = 'two.sided');
}
\arguments{
  \item{x}{Vector of numbers to analyze. Variable must be dichotomous}
  \item{y}{Vector of numbers to analyze}
  \item{alternative}{What is the null-hypothesis?}
}
\value{Returns a two-element vector containing the point-biserial correlation and its p-value.}
\author{Paul C. Boutros & Dorota Sendorek}
\examples{
get.biserial.p.and.corr(
	x = sample(c(0:1), 10, replace= TRUE),
	y = runif(10)
	);
}
\keyword{htest}

\name{get.aov.p}
\alias{get.aov.p}
\title{Get the p-value from a one-way ANOVA}
\description{Returns the p-value from a one-way ANOVA. Returns NA if the t-test is uncalculable (i.e. does try() handling).}
\usage{
get.aov.p(formula, data, ...);
}
\arguments{
  \item{formula}{The ANOVA formula (careful with ordering!)}
  \item{data}{data frame containing data to be fit}
  \item{...}{any extra parameters to be passed to aov()}
}
\value{No return value}
\author{Paul C. Boutros}
\examples{
# need to add an example!
}
\keyword{htest}

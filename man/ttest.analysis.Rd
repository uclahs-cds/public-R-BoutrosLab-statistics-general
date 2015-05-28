\name{ttest.analysis}
\alias{ttest.analysis}
\title{Perform a t-test}
\description{Perform a t-test and return summary statistics}
\usage{
ttest.analysis(values, groups, alternative = "two.sided");
}
\arguments{
  \item{values}{Vector of values to be analyzed}
  \item{groups}{Factor variable with two levels indicating which entries in 'values' belong to which group for the t-test comparison}
  \item{alternative}{Type of alternative hypothesis to be used ("two.sided", "greater", or "less").  Defaults to "two.sided".}
  }
\value{Data frame with relevant statistics}
\author{Nathalie Moon}
\examples{
ttest.analysis(values = c(rnorm(n=50, mean=0), rnorm(n=50, mean=10)), 
    groups = factor(c(rep('A', times=50), rep('B', times=50))));
}
\keyword{htest}

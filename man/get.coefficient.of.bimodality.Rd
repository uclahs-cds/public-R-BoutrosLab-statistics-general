\name{get.coefficient.of.bimodality}
\alias{get.coefficient.of.bimodality}
\title{Coefficient of Biomodality}
\description{Calculates the coefficient of bimodality from a vector, as defined as (1+skewness^2)/(kurtosis+3). For normal distributions, the value should be close to 0.16667; for uniform distributions, close to 0.5555; and for bimodal distributions, the coefficient should be close to 1.}
\usage{
get.coefficient.of.bimodality(x);
}
\arguments{
	\item{x}{A vector of values}
}
\value{A value corresponding to the cofficient of bimodality}
\author{Emilie Lalonde}
\examples{
x <- rnorm(1000);
coeff <- get.coefficient.of.bimodality(x);
}
\keyword{manip}

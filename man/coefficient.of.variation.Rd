\name{coefficient.of.variation}
\alias{coefficient.of.variation}
\title{Calculate coefficient of variation}
\description{Takes a vector and calculates the ratio of the standard deviation to the mean (coefficient of variation, CV)}
\usage{
coefficient.of.variation(x, na.rm);
}
\arguments{
  \item{x}{Vector of numbers to analyze}
  \item{na.rm}{A logical indicating whether or not to ignore (remove) missing values}
}
\value{The modified vector}
\author{Sylvia Shiah}
\examples{
coefficient.of.variation(
	x = rnorm(100)
	);
}
\keyword{manip}

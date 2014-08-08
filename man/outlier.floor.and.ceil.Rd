\name{outlier.floor.and.ceil}
\alias{outlier.floor.and.ceil}
\title{Floor and ceil outliers in a vector}
\description{Takes a vector, identifies a given percentage of outliers on each side, and moves them to the floor or the ceiling respectively.}
\usage{
outlier.floor.and.ceil(x, trim);
}
\arguments{
  \item{x}{Vector of numbers to analyze}
  \item{trim}{The percentage of values to floor/ceil on each side (i.e. total number changed is twice this value)}
}
\value{The modified vector}
\author{Paul C. Boutros}
\examples{
x <- c(-100, rnorm(18), 100);
y <- outlier.floor.and.ceil(x = x, trim = 0.05);
plot(x, y, pch = 19, las = 1);
}
\keyword{manip}

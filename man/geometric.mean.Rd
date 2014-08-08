\name{geometric.mean}
\alias{geometric.mean}
\title{Geometric Mean}
\description{Calculates the Geometric Mean from a vector}
\usage{
geometric.mean(x, logged = FALSE, na.rm = TRUE);
}
\arguments{
  \item{x}{A vector. Values need to be positive if on the log scale.}
  \item{logged}{Indicate whether data is log-transformed or not. Defaults to FALSE. If logged then the function is equivelant to the arithmetic mean.}
  \item{na.rm}{Should NAs be removed.}
}
\value{A value corresponding to the geometric mean}
\author{Daryl Waggott}
\examples{
x <- abs(rnorm(100));
mean1 <- geometric.mean(x);
mean2 <- mean(x, na.rm = TRUE);
median1 <- median(x, na.rm = TRUE);
}
\keyword{manip}

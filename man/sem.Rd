\name{sem}
\alias{sem}
\title{Standard Error of the Mean}
\description{Takes a vector and computes its standard error of the mean.}
\usage{
sem(x, na.rm = TRUE);
}
\arguments{
  \item{x}{A vector}
  \item{na.rm}{Should NAs be removed. Defaults to TRUE.}
}
\value{A numeric value giving the standard error of the mean of the input vector}
\author{Stephenie Prokopec}
\examples{
x <- c(1,2,2,NA,3);
sem(x);
sem(x, na.rm = FALSE);
}
\keyword{manip}

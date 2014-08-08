\name{get.standard.error}
\alias{get.standard.error}
\title{Get the standard error}
\description{Returns the standard error for a vector of numbers. Returns NA if any values in x are NA.}
\usage{
get.standard.error(x);
}
\arguments{
  \item{x}{The numeric vector}
}
\value{The standard error}
\author{Emilie Lalonde}
\examples{
x <- rnorm(100);
get.standard.error(x);
}

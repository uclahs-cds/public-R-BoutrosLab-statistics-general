\name{count.NA}
\alias{count.NA}
\title{Number of missing values in a vector}
\description{Takes a vector and determines the number of missing values in it.}
\usage{
count.NA(x);
}
\arguments{
  \item{x}{A vector}
}
\value{An integer giving the number of NAs in the vector}
\author{Paul C. Boutros}
\examples{
count.NA( c(NA, NA, 1:10) );
}
\keyword{manip}

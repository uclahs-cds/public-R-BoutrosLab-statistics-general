\name{zscore}
\alias{zscore}
\title{zscore transformation or quantile normalization}
\description{Transform a continuous variable into a zscore using ranks.  This is useful for semi-parametric screening in linear models.  It is not suitable for data with many ties}
\usage{
zscore(x);
}
\arguments{
  \item{x}{A vector}
}
\value{A vector of normally distributed values.}
\author{Daryl Waggott}
\examples{
x <- rchisq(100, df = 1);
x.norm <- zscore(x);
}
\keyword{manip}

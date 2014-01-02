\name{get.odds.ratio}
\alias{get.odds.ratio}
\title{Calculate the odds ratio and its confidence}
\description{Takes a pair of binary vectors and computes the log ratio, log-odds, standard error of the log-odds, and 95\% CI of the log-odds.}
\usage{
get.odds.ratio(x,y);
}
\arguments{
  \item{x}{A binary vector}
  \item{y}{A binary vector}
}
\value{A vector containing: odds-ratio, log odds-ratio, standard error of the log odds-ratio, and 95\% CIs of the log-odds ratio}
\author{Paul C. Boutros}
\examples{
x <- sample(0:1, 100, replace = TRUE);
y <- sample(0:1, 100, replace = TRUE);
get.odds.ratio(x,y);
}
\keyword{manip}

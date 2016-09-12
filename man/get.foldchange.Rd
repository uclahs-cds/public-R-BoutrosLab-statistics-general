\name{get.foldchange}
\alias{get.foldchange}
\title{get.foldchange}
\description{Calculates the fold-changes between two samples. Calculated as mean(group2) / mean(group1)}
\usage{
get.foldchange(x, group1, group2, logged = TRUE);
}
\arguments{
  \item{x}{Vector of numbers to analyze}
  \item{group1}{Vector of TRUE/FALSE indicating which are the x samples for the fold-changes}
  \item{group2}{Vector of TRUE/FALSE indicating which are the y samples for the fold-changes}
  \item{logged}{Is the data in log-space or in normal-space?}
}
\value{calculate the fold change of group1 relating to group2}
\author{Sylvia Shiah}
\examples{
get.foldchange(
	x = rnorm(100),
	group1 = c(rep(TRUE, 50), rep(FALSE, 50)),
	group2 = c(rep(FALSE, 50), rep(TRUE, 50))
	);
}

\keyword{htest}


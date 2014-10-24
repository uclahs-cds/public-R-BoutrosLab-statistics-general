\name{get.pve}
\alias{get.pve}
\title{Get the percent variance explained calculated from an analysis of variance table from a fitted model}
\description{Get the percent variance explained (PVE) calculated from an analysis of variance table from a fitted model}
\usage{
get.pve(model);
}
\arguments{
  \item{model}{The model to run the F-test on (e.g. output of lm or glm)}
}
\value{A data-frame with the PVE for each predictor variable (colnames predictor and pve)}
\author{Emilie Lalonde}
\examples{
# need to add an example!
}

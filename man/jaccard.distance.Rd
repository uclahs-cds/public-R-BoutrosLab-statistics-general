\name{jaccard.distance}
\alias{jaccard.distance}
\title{Jaccard distance}
\description{
    Calculates the Jaccard distance between two vectors
    OR Calculates the Jacard distance matrix between several vectors
    The Jaccard distance is the proportion of items that do NOT belong
        to both sets
    }
\usage{
jaccard.distance(x,y = NULL, diag = FALSE, upper = FALSE);
}
\arguments{
  \item{x}{A list of values or binary mask, or matrix}
  \item{y}{A list of values or binary mask or NULL}
  \item{diag}{Whether to include the diagonal of the the distance matrix in output}
  \item{upper}{Whether to include the upper triangular part of the matrix in output}
}
\value{The jaccard distance}
\author{Ryan Kearns}
\examples{
#Binary masks
x <- c(0,1,0,0,0,1,1,1); #1 indicates item belongs to set "x"
y <- c(0,1,0,0,1,0,0,0); #1 indicates item belongs to set "y"
jaccard.distance( x , y );
#Values
x <- c(7.6,3,6,8,3,1,6);
y <- c(6,7,2,6.5,2,0.12,1);
jaccard.distance( x, y );
}
\keyword{manip}

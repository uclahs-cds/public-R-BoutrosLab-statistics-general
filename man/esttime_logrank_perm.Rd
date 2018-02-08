\name{esttime_logrank_perm}
\alias{esttime_logrank_perm}
\title{Estimated run time for Log-Rank Permutation Test}
\description{Calculates the estimated run time for performing a log-rank permutation test with N patients with k variants}
\usage{
esttime_logrank_perm(N,k,no_cores=1,tform='mins');
}
\arguments{
        \item{N}{The number of patients}
        \item{k}{The number of patients in one of the groups}
        \item{no_cores}{Number of cores that will be used to parallelize (efficiency peaks at 6)}
        \item{tform}{The returned time value (default minutes). Options include: secs, mins, hours}
}
\value{Returns the expected run-time. Note, ignores time required to initialize clusters and is only a linear extrapolation. Small errors will magnify for very large number of calculations}
\author{Erik Drysdale}
\examples{
# # Get the estimated time (may take a minute or so but will hold for all N/k)
# N <- 50
# k <- 5
# esttime_logrank_perm(N=N ,k=k, no_cores=4,tform='secs')
# # Check
# So <- matrix(c(1:N,rep(1,N)),ncol=2)
# x <- rep(0,N)
# x[sample(1:N,k)] <- 1
# logrank_perm(So=So, x=x, no_cores=4)
}

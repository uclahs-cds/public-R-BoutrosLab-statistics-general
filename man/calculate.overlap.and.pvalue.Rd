\name{calculate.overlap.and.pvalue}
\alias{calculate.overlap.and.pvalue}
\title{Calculate overlap and p-value}
\description{Calculate actual overlap and corresponding p-value using the hypergeometric distribution as well as the expected overlap.}
\usage{
calculate.overlap.and.pvalue(list1, list2, total.size, lower.tail = TRUE, adjust = FALSE);
}
\arguments{
  \item{list1}{Vector representing a query list}
  \item{list2}{Vector representing a reference list being compared to}
  \item{total.size}{Numerical value representing the total size of reference sample pool}
  \item{lower.tail}{logical for calculating p-value; if TRUE (default), probability is P[X <= x], otherwise, P[X > x] unless adjust = TRUE.  See phyper}
  \item{adjust}{logical for calculations for lower.tail = FALSE; if FALSE (default), P[X > x], otherwise, P[X >= x]}
  }
\details{Calculations are based on using list2 as the basis for comparison but if list1 and list2 come from the same sample pool, it doesn't matter how list1 or list2 are designated. Be careful when interpreting results for lower.tail = FALSE because probability is P[X > x] and not P[X >= x]. }
\value{Returns a three-element vector containing the actual overlap, expected overlap and p-value}
\author{Paul C. Boutros and Denise Mak}
\examples{
# SNP array example 1
# Total number of SNPs on array is 2000
# Identified 56 special SNPs in group1
# Identified 78 special SNPs in group2
# 10 interesting SNPs are common to both groups
# Question: Is 10 greater or less than what we would expect?

list1 <- paste("SNP", 1:56);
list2 <- paste("SNP", 47:124);
total.size <- 2000;

# Expected overlap is 2.2

# SNP array example 2
# Total number of SNPs on array is 2000
# 500 SNPs on array found on chromosome 1 
# Identified 41 special SNPs
# 5 interesting SNPs found on chromosome 1
# Question: Is 5 over or under-represented on chromosome 1 as compared to
# the proportion of chromosome 1 SNPs on the array?

list1 <- paste("SNP", 1:41);
list2 <- paste("SNP", 37:536);
total.size <- 2000;

# Expected overlap is 10.3

# Calculate if over-represented
calculate.overlap.and.pvalue(
	list1 = list1, 
	list2 = list2,
	total.size = total.size,
	lower.tail = FALSE,
	adjust = TRUE
	);
# Note: calculation is for P[X > x]

# Calculate if under-represented
calculate.overlap.and.pvalue(
	list1 = list1, 
	list2 = list2,
	total.size = total.size,
	lower.tail = TRUE,
	adjust = FALSE
	);
# Note: calculation is for P[X <= x]
}
\seealso{
  \code{\link{phyper}}
  }
\keyword{htest}

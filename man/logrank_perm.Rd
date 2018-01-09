\name{logrank_perm}
\alias{logrank_perm}
\title{Log-Rank Permutation Test}
\description{Calculates log-rank statistic for all (N k) possible outcomes and associated p-value under null of no survival differences}
\usage{
logrank_perm(So,x,no_cores=1);
}
\arguments{
        \item{So}{Matrix-formatted Surv() object}
        \item{x}{Binary vector for two groups}
        \item{no_cores}{Number of cores to use}
}
\value{A value corresponding to the cofficient of bimodality}
\author{Erik Drysdale}

\examples{
# # WARNING! FOR LARGE NUMBERS OF CALCULATIONS, REQUEST SUFFICIENT RAM OF PROGRAM MAY CRASH!
# # Matrix version
# So <- as.matrix(survival::Surv(time=survival::aml$time,event=survival::aml$status))
# # Encode one of the groups
# x <- ifelse(survival::aml$x=='Maintained',1,0)
# # Run the classical log-rank test
# lrtest <- survival::survdiff(survival::Surv(time=survival::aml$time,event=survival::aml$status)~x)
# # Permutation version (may take a minute)
# permtest <- logrank_perm(So=So,x=x,no_cores = 12)
# # Compare p-values
# data.frame(lrtest=pchisq(lrtest$chisq,1,lower.tail=FALSE),
# 		   permtest=permtest$pval)
}

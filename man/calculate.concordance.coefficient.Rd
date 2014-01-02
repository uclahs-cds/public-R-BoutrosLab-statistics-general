\name{calculate.concordance.coefficient}
\alias{calculate.concordance.coefficient}
\title{Calculate Lin's concordance correlation coefficient}
\description{Calculates the concordance correlation coefficient as derived by Lin for a set of matched observations in multiple classes.}
\usage{
calculate.concordance.coefficient(data.mat);
}
\arguments{
  \item{data.mat}{A matrix of numeric data, with columns representing classes and row representing matched observations}
}
\details{This function calculates the Lin concordance correlation coefficient for a set of matched observations for any number of observers. Roughly speaking, this metric measures the orthogonal distance of the data from the y=x line. For further information on Lin's methods, see 'A concordance correlation coefficient to evaluate reproducability', Biometrics 1989, 45:1, pp. 255-68 (Section 7) and 'A note on the concordance correlation coefficient', Biometrics 2000, 56:1, pp. 324-325 (A correction to the formula given in the former).}
\value{Returns a single value on the interval [-1,1] representing the degree to which the classes are in concordance with one another (interpreted the same as a Pearson or a Spearman correlation coefficient).}
\author{Michal Grzadkowski}
\examples{
# a matrix with four classes of 100 random observations each
data.matrix1 <- matrix(
			runif(400),
			ncol = 4
			);
calculate.concordance.coefficient(data.matrix1);

# a matrix with four classes of 100 similar observations each
base.data <- sample(
		x = 0:1000,
		size = 100,
		replace = TRUE
		);
data.matrix2 <- matrix(
			c(
				base.data,
				base.data + runif(100),
				base.data + runif(100),
				base.data + runif(100)
				),
			ncol = 4
			);
calculate.concordance.coefficient(data.matrix2);
}


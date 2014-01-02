\name{weighted.AIC}
\alias{weighted.AIC}
\title{Weighted AIC}
\description{Calculate weighted Akaike Information Criterion scores}
\usage{
weighted.AIC(x);
}
\arguments{
	\item{x}{A vector of AIC values, must have greater than 2 elements}
}
\value{A vector of weighted AIC values corresponding to the input. }
\details{Weighted AIC values can be interpreted as the probability of each model being the best given the data and thet set of models being compared. }
\references{More information including formula for weighted AIC can be found in: Eric-Jan Wagenmakers and Simon Farrell. AIC model selection using Akaike weights. Psychonomic Bulletin and Review, 11(1):192-6, Feb 2004.}
\author{Denise Mak}
\examples{
set.seed(123456789);

# creating sample non-linear data
temp <- data.frame(x = 1:100);
temp$y <- temp$x^2-rnorm(100, 0, 500);

# fitting with different models
# linear model
fit1 <- lm(y~x, data = temp);
# quadratic model
fit2 <- lm(y~x + I(x^2), data = temp);
# logistic model
fit3 <- nls(y~SSlogis(x, Asym, xmid, scal), data = temp);

\dontrun{
# Example to visualize different fits
plot(temp);
lines(temp$x, predict(fit1), col = "blue")
lines(temp$x, predict(fit2), col = "red")
lines(temp$x, predict(fit3), col = "green")
}

# calculate AIC values (see stats::AIC for more information)
x <- AIC(fit1, fit2, fit3);

# quadratic model has the highest probability and the best fit
weighted.AIC(x$AIC);
}
\seealso{
  \code{\link{AIC}}
  }
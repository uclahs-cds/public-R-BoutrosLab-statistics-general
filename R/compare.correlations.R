# The BoutrosLab.statistics.general package is copyright (c) 2011 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

compare.correlations <- function(corr1, corr2, n1, n2, method) {

	if ('spearman' == method) {

		# apply Fisher transform on the correlation coefficients
		rho1.transformed  <- 0.5 * log( (1+corr1) / (1-corr1) );
		rho2.transformed  <- 0.5 * log( (1+corr2) / (1-corr2) );
		delta.rho         <- rho2.transformed - rho1.transformed;
	
		# calculate the standard deviation of spearman correlation coefficients
		# the formula and the value of parameter **a** is due to Fieller et al.
		# "Tests for Rank Correlation Coefficients. I," Biometrica, Vol. 44, No. 3/4 (Dec. 1957)
		a    <- 1.060;
		sd1  <- a / (n1 - 3);
		sd2  <- a / (n2 - 3);
		sd12 <- sqrt(sd1 + sd2);

		p.value <- 2 * pnorm(
			abs(delta.rho),
			mean = 0,
			sd = sd12,
			lower.tail = FALSE
			);

		# return as a vector to allow easier usage in apply() contexts
		return(	c(corr1, corr2, sd1, sd2, delta.rho, sd12, p.value, n1, n2) );

		}

	else {
		print('Invalid correlation method');
		return(NA);
		}

	}

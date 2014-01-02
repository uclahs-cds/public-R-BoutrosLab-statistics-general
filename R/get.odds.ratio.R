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

get.odds.ratio <- function(x,y) {

	# verify that the length of the two variables is equivalent
	if (length(x) != length(y)) {
		warning('Non-matching lengths');
		return(NA);
		}

	# verify that we have binary variables
	if (2 != nlevels(as.factor(x)) & 2) {
		warning('Non binary x variable');
		return(NA);
		}
	else if (2 != nlevels(as.factor(y))) {
		warning('Non binary y variable');
		return(NA);
		}

	# construct a two-way table
	z <- table(x,y);

	# calculate the odds ratio from the joint distribution
	odds.ratio <- z[1,1] * z[2,2] / z[1,2] / z[2,1];

	# calculate the standard error
	standard.error <- sqrt(1/z[1,1] + 1/z[1,2] + 1/z[2,1] + 1/z[2,2]);

	# get the 95% CIs
	U95 <- log(odds.ratio) + 1.96 * standard.error;
	L95 <- log(odds.ratio) - 1.96 * standard.error;

	# return stats to the caller
	return(
		c(
			odds.ratio,
			log(odds.ratio),
			standard.error,
			L95,
			U95
			)
		);

	}

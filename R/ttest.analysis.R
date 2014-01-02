# The BoutrosLab.statistics.general package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

ttest.analysis <- function(values, groups, alternative = "two.sided"){

	# Check that 'groups' only takes on two values, as a t-test is limited to the comparison
	# of two groups
	if (length(unique(groups)) > 2){
		stop('ttest is limited to comparison of two groups, and you specified more than two groups');
		}

	# Ensure that 'groups' is a factor, so that groups are refered to in the appropriate order
	groups <- as.factor(groups);

	# perform t-test to assess mean difference between groups
	stats <- t.test(formula = values ~ groups, alternative = alternative);

	# define statistical parameters
	this.pval <- stats$p.value;
	average.group1 <- stats$estimate[1];
	average.group2 <- stats$estimate[2];

	# prepare summary statistical values to return, if requested
	summary.stats <- data.frame(
		statistical.method = "ttest",
		alternative.hypothesis = alternative,
		pvalue = this.pval,
		group.1 = levels(groups)[1],
		group.2 = levels(groups)[2],
		mean.group.1 = average.group1,
		mean.group.2 = average.group2
		);

	# Remove ret.stats rownames to make returned value look nicer
	rownames(summary.stats) <- NULL;

	# Return summary.stats
	return(summary.stats)
	}

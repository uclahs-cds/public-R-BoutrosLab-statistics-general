# The BoutrosLab.statistics.general package is copyright (c) 2016 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.


### GET.TOP.OUTLIERS ##########################################################################
# Compute rho of 2 datasets, with or without outliers.
# Input variable:
#  - data frame of 2 datasets organized in columns with meaningful row names (e.g. gene names)
#  - max.outliers, number of outliers to skip
#  - optional input: delta.limit, minimum change in correlation coefficient for a data point to be considered an outlier. Default to 0.05
#  - optional input: n, number of outlier combinations to report
#  - default value n = 5, report 5 improvement
#  - optional input: correlation method. Default to spearman.
# Output variable:
#  - outputs list of top n new correlation coefficients in descending order, with outlier samples by row name, along with original coefficient in first row


get.top.outliers <- function(input.df, max.outliers, delta.limit = 0.05, n = 5, cor.type = 'spearman') {
	if (!is.data.frame(input.df) | !is.numeric(max.outliers) | !is.numeric(delta.limit) | !is.numeric(n)) {
		# check if input format is valid
		stop('Invalid input format.');
	} else if (nrow(input.df) < max.outliers | max.outliers < 0) {
		# check if max.outliers has a valid value
		stop('Invalid outlier number. Make sure it\'s less than total sample number.');
	} else if ('pearson' != cor.type & 'kendall' != cor.type & 'spearman' != cor.type) {
		# check correlation method is valid
		stop('Invalid correlation method. Accepted methods: \'pearson\', \'kendall\', \'spearman\'.');
	} else {
		corr.unadj <- unname(
			cor.test(
				input.df[,1], 
				input.df[,2], 
				method = cor.type
				)$estimate
			);
		# screen data, choose only data points that will cause a change in rho (or any other correlation value)
		corr.prelim <- apply(
			matrix(seq(1, nrow(input.df), 1)),
			1,
			function(x) unname(
				cor.test(
					input.df[-x,1], 
					input.df[-x,2],
					method = cor.type
					)$estimate
				)
			);
		signi.samples <- which(corr.prelim > (corr.unadj + delta.limit));		
		# permutation of all possible combinations of 'n' samples to skip within list of significant samples
		# first make sure that n is less than or equal to number of significant samples
		if (0 == length(signi.samples)) {
			output.sample <- data.frame(
				c(corr.unadj, 0), 
				c(0, 0), 
				replicate(max.outliers, c(NA, NA)));
			colnames(output.sample) <- c('new.corr', 'corr.delta', paste('outlier.sample', 1:max.outliers, sep = '.'));
			return(output.sample);
		} else if (max.outliers > length(signi.samples)) {
			max.outliers = length(signi.samples);
			};
		outlier.samples <- data.frame(
			matrix(
				combn(signi.samples, max.outliers), 
				ncol= max.outliers,
				byrow = TRUE
				)
			);
		if (n > nrow(outlier.samples)) {
			n = nrow(outlier.samples)
		};
		# compute correlation coefficient for each combination of outliers
		corr.adj <- apply(
			outlier.samples,
			1,
			function(x) unname(
				cor.test(
					input.df[-x,1], 
					input.df[-x,2],
					method = cor.type
					)$estimate
				)
			);
		corr.delta <- corr.adj - corr.unadj;
		output.sample <- data.frame(corr.adj, corr.delta, matrix(rownames(input.df)[unlist(outlier.samples)], ncol=max.outliers));
		colnames(output.sample) <- c('new.corr', 'corr.delta', paste('outlier.sample', 1:max.outliers, sep = '.'));
		# sort by descending corr delta
		output.sample <- output.sample[order(-output.sample$corr.delta),];
		output.sample <- rbind(c(corr.unadj, 0, rep(NA, max.outliers)), output.sample[c(1:n), ])
		};
	return(output.sample);
	};
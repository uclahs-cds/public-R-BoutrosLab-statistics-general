# The BoutrosLab.statistics.general package is copyright (c) 2010 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

calculate.and.compare.correlations <- function(x1, x2, y1, y2, method = 'spearman') {

	# ensure a viable correlation method is specified
	if (! method %in% c('spearman')) {
		print('Invalid correlation method');
		return(NA);
		}

	# calculate the correlations
	corr1 <- cor(x1, y1, method = method, use = "pairwise.complete.obs");
	corr2 <- cor(x2, y2, method = method, use = "pairwise.complete.obs");

	# calculate the significance of the correlation
	comparison <- compare.correlations(
		corr1 = corr1,
		corr2 = corr2,
		n1 = length(x1[!is.na(x1) & !is.na(y1)]), #calculating jointly non-NA length 
		n2 = length(x2[!is.na(x2) & !is.na(y2)]), #calculating jointly non-NA length 
		method = method
		);

	return(comparison);

	}

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

get.ttest.p.and.foldchange <- function(x, group1, group2, paired = FALSE, var.equal = FALSE, alternative = 'two.sided', logged = TRUE) {
	
	group1 <- as.numeric(x[group1]);
	group2 <- as.numeric(x[group2]);

	tryCatch(
		expr = c(
			t.test(
				group1,
				group2,
				paired = paired,
				var.equal = var.equal,
				alternative = alternative
				)$p.value,
			if (logged == TRUE) { mean(group2, na.rm = TRUE) - mean(group1, na.rm = TRUE); } 
			else { mean(group2, na.rm = TRUE) / mean(group1, na.rm = TRUE); }
			),
		 error = function(e) { return ( c(NA, NA)); }
		 );
	}

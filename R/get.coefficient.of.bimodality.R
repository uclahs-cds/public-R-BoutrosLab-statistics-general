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

get.coefficient.of.bimodality <- function (x){
	
	# get mean and standard deviatino of data
	meanX <- mean(x);
	sdev <- sd(x);
	
	# explicitly give dimensions to data
	dataX <- array(x, dim=length(x));	
	
	
	# calculate the skewness
	dif <- apply (dataX, 1, function(a) (a-meanX)^3);
	skewness <- sum(dif) / ( (length(x)-1) * (sdev^3) )	
	
	# calculate the kurtosis
	dif2 <- apply (dataX, 1, function(a) (a-meanX)^4);
	kurtosis <- sum(dif2) / ( (length(x)-1) * (sdev^4) );	
	
	# calculate the coefficient of bimodality
	coeff.of.bimod <- (1 + (skewness^2) ) / (kurtosis +3);
	return(coeff.of.bimod);
	
	}

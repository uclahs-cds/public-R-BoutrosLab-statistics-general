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

calculate.concordance.coefficient <- function(data.mat) {

        # gets the covariance matrix and the means of each column of observations
        covariance.mat <- cov(data.mat);
        mean.vec <- colMeans(data.mat);
        mean.mat <- outer(mean.vec, mean.vec, '-') ^ 2;

        # calculates Lin's CCC and returns it
        ccc.value <- (2 * sum(covariance.mat[upper.tri(covariance.mat)])) / ((sum(mean.mat[upper.tri(mean.mat)]) / 2) + ((ncol(data.mat) - 1) * sum(diag(covariance.mat))));
        return(ccc.value);

        }


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

jaccard.distance <- function(x, y = NULL, diag = FALSE, upper = FALSE) {

	if (is.data.frame(y)) { y <- as.matrix(y); }
    if (is.data.frame(x)) { x <- as.matrix(x); }

    if (!is.matrix(x) && is.null(y)) {
        stop("supply both 'x' and 'y' or a matrix-like 'x'");
		}
	else if (!(is.numeric(x) || is.logical(x))) {
        stop("'x' must be numeric");
		}

	stopifnot(is.atomic(x));

    if (!is.null(y)) {

		if (!(is.numeric(y) || is.logical(y))) { stop("'y' must be numeric"); }
        stopifnot(is.atomic(y));

        if (dim(as.matrix(x))[1] != dim(as.matrix(y))[1]) {
            stop("'x' and 'y' must be same length");
			}

		x <- cbind(x,y);
		}

    if (requireNamespace("vegan")) {
		return(vegan::vegdist(x,diag = diag, upper = upper, method = 'jaccard')); #use C version
		}

    datalists <- dim(x)[1];
    length <- dim(x)[2];

    if (100 < datalists) {
		cat("Warning: Without vegan package, jaccard.distance runs very slowly on larger datasets.\neg. A 1000x1000 dataset takes 30 minutes");
		}

    outmatrix <- mat.or.vec(datalists,datalists);

    for (j in 1:(datalists-1)) {
        for (i in (j+1):datalists) {
			A <- sum(x[i,]) + sum(x[j,]);
            J <- sum( apply(cbind(x[i,], x[j,]),1, min) );
            outmatrix[i,j] <- (A-2*J)/(A-J);
			}
		}

    if (upper){
		outmatrix <- outmatrix + t(outmatrix); #fill in upper portion
		}

	outmatrix <- as.dist(outmatrix, diag = diag, upper = upper)

	return(outmatrix);

	}


get.pve  <- function(model){
	anova.out <- anova(model);

	pve <- anova.out$Deviance[-1] / sum(anova.out$Deviance[-1], na.rm = TRUE);
	out.df <- data.frame(
		predictor = rownames(anova.out),
		pve = pve
		);
	return(out.df);
	}
